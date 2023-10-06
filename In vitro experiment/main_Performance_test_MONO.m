% Imatinib-sensitive and -resistant Ba/F3 cells DATA
    load('MONOCLONAL_DATA.mat')
    DATA = SENSITIVE_500_BF; % Selecting the data from [BF_11,BF_12,BF_21,BF_41].
    init = 500;
    [NR,NC,NT] = size(DATA);
    lower_mixture_limit = 0.10;

    Time = [9:3:48];
    Conc = 10^(6)*[0  31.25*10^(-9) 62.5*10^(-9) 125*10^(-9) 250*10^(-9) 375*10^(-9) 500*10^(-9) 1.25*10^(-6) 2.5*10^(-6) 3.75*10^(-6) 5*10^(-6) ];
%     Conc = 10^(6)*[0  31.25*10^(-9) 62.5*10^(-9) 125*10^(-9) 250*10^(-9) 500*10^(-9) 750*10^(-9) 1*10^(-6) 2.5*10^(-6) 5*10^(-6) 7.5*10^(-6) ];

    alpha_ub = 0.06;
    alpha_lb = 1e-6;
    b_ub = 1;
    b_lb = 0.878;
    E_ub = 50;
    E_lb = 1e-6;
    n_ub = 20;
    n_lb = 1e-3;
    mix_ub = 1;
    sig_ub = 2500;

    % BD initial points
    lb_BD=zeros(5,1);
    ub_BD=lb_BD;
    ub_BD(1)=1; % natural birth rate
    ub_BD(2)=1; % natural death rate
    ub_BD(3)=1; % b
    lb_BD(3)=0.878;
    ub_BD(4)=50; % E
    lb_BD(4)=1e-6;
    ub_BD(5)=20; % n
    lb_BD(5)=1e-3;
    lb_BD   = [lb_BD;0];
    ub_BD   = [ub_BD;100];
    A_BD    = [0,1,-1,0,0,0;
               0,-1,1,0,0,0];
    b_BD    = [0.06;0];
    
    num_sub = 1;
    lb_HL= zeros(5*num_sub + 1,1);
    ub_HL=sig_ub*ones(size(lb_HL));
    
    ub_HL(1:5:5*num_sub-4)=alpha_ub; % alpha
    ub_HL(2:5:5*num_sub-3)=b_ub;  % b
    ub_HL(3:5:5*num_sub-2)=E_ub; % E
    ub_HL(4:5:5*num_sub-1)=n_ub; % n
    ub_HL(5:5:5*num_sub-5)=mix_ub; % Mixture parameter
    
    lb_HL(1:5:5*num_sub-4)=alpha_lb; % alpha
    lb_HL(2:5:5*num_sub-3)=b_lb;  % b
    lb_HL(3:5:5*num_sub-2)=E_lb; % E
    lb_HL(4:5:5*num_sub-1)=n_lb; % n

    A_HL=zeros(1,length(ub_HL));
    A_HL(5:5:5*num_sub-5) = 1;
    b_HL = 1 - lower_mixture_limit;
    
    ConcT = 1;      
    TimeT = 30; 





    %% Optimization
    

    num_optim = 100;
    options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',500,'Display','iter','algorithm','sqp');
%     func3  = @(x) Cov_diff(DATA, Conc, Time, x,2); 
    func_hl   = @(x)HL_mono(x,init,DATA,Conc,Time, ConcT, TimeT, NR);
    x_init_HL    = [];
    for n=1:num_optim
        x02      = rand(length(ub_HL),1).*(ub_HL-lb_HL) + lb_HL;
    %     x02      = [x02(1);x02(3:12)];
        x_init_HL   = [x_init_HL,x02];
    end
    
    

    %% High Low variance model
    fval2_hl   = [];
    params2_hl = [];
    grad2_hl   = [];
    time2_hl   = [];
    
    parfor n=1:num_optim
        [xx_hl,ff_hl,~,out_hl,~,g_hl,~]  = fmincon(func_hl,x_init_HL(:,n),A_HL,b_HL,[],[],lb_HL,ub_HL,[],options1);
    %     grad     = norm(g,inf);
        fval2_hl    = [fval2_hl, ff_hl];
        params2_hl  = [params2_hl,xx_hl];
        grad2_hl    = [grad2_hl,g_hl];
    end
%     [opt_fval_hl, opt_xx_hl] = vectorized_inference_k_subpopulations_two_noise_levels(DATA, Conc, ConcT, TimeT, ub_HL, lb_HL, num_optim, num_sub, A_HL, b_HL);



    %% BD model
    x_init_BD    = [];

    for n=1:num_optim
        x02      = rand(length(ub_BD),1).*(ub_BD-lb_BD) + lb_BD;
    %     x02      = [x02(1);x02(3:12)];
        x02(2) = max(0,x02(1) - rand*0.1);
        x_init_BD   = [x_init_BD,x02];
    end


    
    %% Dynamic variance model


%%  dyn optimization
    func_dyn  = @(x) BD_Obj_mono(x,DATA,Conc,Time,NR,init);
    func_sto  = @(x) BD_obj_func_mono(DATA, Conc, Time, NR, x,init);
    fval2_dyn   = [];
    params2_dyn = [];
    grad2_dyn   = [];
    
    parfor n=1:num_optim
        [xx_dyn,ff_dyn,~,out_dyn,~,g_dyn,~]  = fmincon(func_dyn,x_init_BD(:,n),A_BD,b_BD,[],[],lb_BD,ub_BD,[],options1);
    %     grad     = norm(g,inf);
        fval2_dyn    = [fval2_dyn, ff_dyn];
        params2_dyn  = [params2_dyn,xx_dyn];
        grad2_dyn    = [grad2_dyn,g_dyn];
    end



    %% Approximated Stochastic model
    fval2_sto  = [];
    params2_sto = [];
    grad2_sto = [];
    
    parfor n=1:num_optim
        [xx_sto,ff_sto,~,out_sto,~,g_sto,~]  = fmincon(func_sto,x_init_BD(:,n),A_BD,b_BD,[],[],lb_BD,ub_BD,[],options1);
        grad     = norm(g_sto,inf);
        fval2_sto    = [fval2_sto, ff_sto];
        params2_sto  = [params2_sto,xx_sto];
        grad2_sto    = [grad2_sto,g_sto];
    end

    %% Get optimal
    [opt_fval_hl, idx_hl]   = min(fval2_hl);
    [opt_fval_dyn, idx_dyn] = min(fval2_dyn);
    [opt_fval_sto, idx_sto] = min(fval2_sto);
    opt_xx_dyn = params2_dyn(:,idx_dyn);
    opt_xx_sto = params2_sto(:,idx_sto);
    opt_xx_hl  = params2_hl(:,idx_hl);
    
    
    
    %% Compute indi
    sto_indi = get_indi_mono(opt_xx_sto,max(Conc),0);
    dyn_indi = get_indi_mono(opt_xx_dyn,max(Conc),0);
    hl_indi  = get_indi_mono(opt_xx_hl,max(Conc),1);
    
    
    %% Compute the AIC/BIC
    
    aic_hl  = compute_aic(opt_fval_hl,length(opt_xx_hl));
    aic_dyn = compute_aic(opt_fval_dyn,length(opt_xx_dyn));
    aic_sto = compute_aic(opt_fval_sto,length(opt_xx_sto));


    save('Result\SENSITIVE_500_HL_MONO.mat')