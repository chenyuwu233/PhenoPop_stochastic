%% Generation of 3 sub-population
rng(2)

%% Initialization
    
    options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',500,'Display','off','algorithm','sqp');
    num_optim= 100;
    num_sub = 3;
    lower_mixture_limit = 0.10;
    initial = 1000;
    NT      = 13; 
    NC      = 11;
    NR      = 20;   
    Conc = 10^(6)*[0 31.25*10^(-9) 62.5*10^(-9) 125*10^(-9) 250*10^(-9) 375*10^(-9) 500*10^(-9) 1.25*10^(-6) 2.5*10^(-6) 3.75*10^(-6) 5*10^(-6)];
    Time =   3*(0:NT-1);

    %% Set feasible region
    beta_ub  = 1;
    beta_lb  = 0;
    nu_ub  = 1;
    nu_lb  = 0;
    alpha_ub = 0.1;
    alpha_lb = 1e-6;
    b_ub = 1-1e-6;
    b_lb = 0.27;
    E_ub = 10;
    E_lb = 1e-6;
    n_ub = 10;
    n_lb = 1e-3;
    mix_ub = 1;
    on_lb = 1;
    on_ub = 10;
    sig_ub = 2500;


    lb_BD = zeros(6*num_sub,1);
    ub_BD = lb_BD;
    lb_BD(end) = on_lb;
    ub_BD(end) = on_ub;

    A_BD = zeros(2*num_sub,6*num_sub);
    b_BD = zeros(2*num_sub,1);
    for i = 1:num_sub
        % theta
        if i == num_sub
            ub_BD(6*i-5) = beta_ub;
            ub_BD(6*i-4) = nu_ub;
            ub_BD(6*i-3) = b_ub;
            ub_BD(6*i-2) = E_ub;
            ub_BD(6*i-1) = n_ub;

            lb_BD(6*i-5) = beta_lb;
            lb_BD(6*i-4) = nu_lb;
            lb_BD(6*i-3) = b_lb;
            lb_BD(6*i-2) = E_lb;
            lb_BD(6*i-1) = n_lb;
            
            % A,b
            beta_nu       = zeros(1,6*num_sub);
            nu_beta       = zeros(1,6*num_sub);
            beta_nu(i*6-5:i*6-4) = [1,-1];
            nu_beta(i*6-5:i*6-4) = [-1,1];
            A_BD(2*i-1,:) = beta_nu;
            A_BD(2*i,:)   = nu_beta;
            b_BD(2*i-1)   = alpha_ub;
            b_BD(2*i)     = 0;
        else
            ub_BD(6*i-5) = mix_ub;
            ub_BD(6*i-4) = beta_ub;
            ub_BD(6*i-3) = nu_ub;
            ub_BD(6*i-2) = b_ub;
            ub_BD(6*i-1) = E_ub;
            ub_BD(6*i)   = n_ub;

            lb_BD(6*i-5) = 0;
            lb_BD(6*i-4) = beta_lb;
            lb_BD(6*i-3) = nu_lb;
            lb_BD(6*i-2) = b_lb;
            lb_BD(6*i-1) = E_lb;
            lb_BD(6*i)   = n_lb;

            % A,b
            beta_nu       = zeros(1,6*num_sub);
            nu_beta       = zeros(1,6*num_sub);
            beta_nu(i*6-4:i*6-3) = [1,-1];
            nu_beta(i*6-4:i*6-3) = [-1,1];
            A_BD(2*i-1,:) = beta_nu;
            A_BD(2*i,:)   = nu_beta;
            b_BD(2*i-1)   = alpha_ub;
            b_BD(2*i)     = 0;
        end


    end
    
    if num_sub >2
        ip_const = zeros(1,6*num_sub);
        ip_const(1:6:6*num_sub-6) = 1;

        A_BD = [A_BD;ip_const];
        b_BD = [b_BD;0.9];
    end




    
    % HL initial
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
    TimeT = 18; 


    x_init_HL    = [];
    for n=1:num_optim
        x02      = rand(length(ub_HL),1).*(ub_HL-lb_HL) + lb_HL;
        if num_sub > 2
            for s = 2:num_sub
                x02(5*s-5) = rand*(1- sum(x02(5:5:5*(s-1)-5)));
            end
        end
        x_init_HL   = [x_init_HL,x02];
    end
    x_init_BD    = [];
    
    for n=1:num_optim
        x02      = rand(length(ub_BD),1).*(ub_BD-lb_BD) + lb_BD;
    
        for s = 1:num_sub
            if s == num_sub
                x02(6*s-4) = max(0,x02(6*s-5)-rand*alpha_ub);
            else
                x02(6*s-3) = max(0,x02(6*s-4)-rand*alpha_ub);
            end
        end
        if num_sub > 2
            for s = 1:num_sub-2
                x02(1+6*s) = rand*(1-sum(x02(1:6:6*s-5)));
            end
        end
        x_init_BD   = [x_init_BD,x02];
    end


    %% Generation parameters
    
    num_sub_GE = 3;
    
    lb_GE = lb_BD(6*(num_sub-num_sub_GE)+1:end);
    ub_GE = ub_BD(6*(num_sub-num_sub_GE)+1:end);
    lb_GE(1:6:end-6) = 1/(2*num_sub_GE);
    ub_GE(1:6:end-6) = 1/(num_sub_GE);
    lb_GE(4:6:end-3) = 0.8;
    ub_GE(4:6:end-3) = 0.9;
    lb_GE(6:6:end-1) = 1.5;
    ub_GE(6:6:end-1) = 5;
    lb_GE(end-3)     = 0.8;
    ub_GE(end-3)     = 0.9;
    lb_GE(end-1)     = 1.5;
    ub_GE(end-1)     = 5;
    
    conc_int = floor(length(Conc)/num_sub_GE);
    for i = 1:num_sub_GE
        E_GE_lb = Conc(2+conc_int*(i-1));
        E_GE_ub = Conc(3+conc_int*(i-1));
        if i~=num_sub_GE
            lb_GE(6*i-1) = E_GE_lb;
            ub_GE(6*i-1) = E_GE_ub;
        else
            lb_GE(6*i-2) = E_GE_lb;
            ub_GE(6*i-2) = E_GE_ub;
        end
    end
    
    
    
    theta = rand(length(ub_GE),1).*(ub_GE-lb_GE) + lb_GE;
    theta(3:6:end-4) = max(0,theta(2:6:end-5) - rand(num_sub_GE-1,1) * 0.1);
    theta(end-4) = max(0,theta(end-5)-rand*0.1);
    
    indi_ip = GR_sort_multi(theta,Conc(end),num_sub_GE);
    
    DATA = sto_gen_bd(NR,Conc,Time,initial,theta,num_sub_GE);



%% Performing Bootstrapping

    B_num = 13;
    B_sample =100;
    Boot_hl  = [];
    Boot_dyn = [];
    Boot_sto = [];
    Boot_hl_indi  = [];
    Boot_dyn_indi = [];
    Boot_sto_indi = [];
    B_indi   = [];

    parfor i = 1:B_sample
                
                Data_i =[];
                indi = [];
                for j = 1:B_num
                    ri = randi(20);
                    indi = [indi,ri];
                    Data_i = [Data_i;DATA(ri,:,:)];
                end
                B_indi    = [B_indi;indi];
            
                func_dyn  = @(x) BD_Obj(x,Data_i,Conc,Time,B_num,num_sub_GE);
                func_hl   = @(x)vectorized_objective_function_two_noise_levels_k_subpopulations(x,num_sub_GE,Data_i,Conc,Time, ConcT, TimeT, B_num);
                func_sto  = @(x) BD_obj_func(Data_i, Conc, Time, B_num, x,num_sub_GE);
    
    
                %% High Low variance model
                fval2_hl   = [];
                params2_hl = [];
                grad2_hl   = [];
                time2_hl   = [];
                
                for n=1:num_optim
                    [xx_hl,ff_hl,~,out_hl,~,g_hl,~]  = fmincon(func_hl,x_init_HL(:,n),A_HL,b_HL,[],[],lb_HL,ub_HL,[],options1);
                %     grad     = norm(g,inf);
                    fval2_hl    = [fval2_hl, ff_hl];
                    params2_hl  = [params2_hl,xx_hl];
                    grad2_hl    = [grad2_hl,g_hl];
                end
        
        
                %% Dynamic variance model
                fval2_dyn   = [];
                params2_dyn = [];
                grad2_dyn   = [];
                
                for n=1:num_optim
                    [xx_dyn,ff_dyn,~,out_dyn,~,g_dyn,~]  = fmincon(func_dyn,x_init_BD(:,n),A_BD,b_BD,[],[],lb_BD,ub_BD,[],options1);
                %     grad     = norm(g,inf);
                    fval2_dyn    = [fval2_dyn, ff_dyn];
                    params2_dyn  = [params2_dyn,xx_dyn];
                    grad2_dyn    = [grad2_dyn,g_dyn];
                end
        
        
                %% Stochastic Model
        
                fval2_sto  = [];
                params2_sto = [];
                grad2_sto = [];
                
                for n=1:num_optim
                    [xx_sto,ff_sto,~,out_sto,~,g_sto,~]  = fmincon(func_sto,x_init_BD(:,n),A_BD,b_BD,[],[],lb_BD,ub_BD,[],options1);
                    grad     = norm(g_sto,inf);
                    fval2_sto    = [fval2_sto, ff_sto];
                    params2_sto  = [params2_sto,xx_sto];
                    grad2_sto    = [grad2_sto,g_sto];
                end
        
        
                %% Obtain the optimization inference
                    
                [opt_fval_hl, idx_hl]   = min(fval2_hl);
                [opt_fval_dyn, idx_dyn] = min(fval2_dyn);
                [opt_fval_sto, idx_sto] = min(fval2_sto);
                opt_xx_hl  = params2_hl(:,idx_hl);
                opt_xx_dyn = params2_dyn(:,idx_dyn);
                opt_xx_sto = params2_sto(:,idx_sto);
                Boot_hl  = [Boot_hl,opt_xx_hl];
                Boot_dyn = [Boot_dyn,opt_xx_dyn];
                Boot_sto = [Boot_sto,opt_xx_sto];
                Indi_xx_hl  = GR_sort_multi(opt_xx_hl,Conc(end),num_sub_GE);
                Indi_xx_dyn = GR_sort_multi(opt_xx_dyn,Conc(end),num_sub_GE);
                Indi_xx_sto = GR_sort_multi(opt_xx_sto,Conc(end),num_sub_GE);
                Boot_hl_indi  = [Boot_hl_indi,Indi_xx_hl'];
                Boot_dyn_indi = [Boot_dyn_indi,Indi_xx_dyn'];
                Boot_sto_indi = [Boot_sto_indi,Indi_xx_sto'];
    
    
    end


save('Result\CI_multi.mat')

