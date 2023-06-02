%% This script is for computing the confidence interval of different parameter generated from the different model.




%% Set Parallel 

parpool('local',24)


%% Setting the seed
rng(111)

% load("Confidence_Interval40(New_gen_dn_sample_long_time).mat")
%% Set up the feasible region
    
    options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',500,'Display','iter','algorithm','sqp');
    num_optim= 100;
    num_sub = 2;
    lower_mixture_limit = 0.10;
    initial = 1000;
    %% HL initial points
        alpha_ub = 0.10;
        alpha_lb = 1e-6;
        b_ub = 1-1e-6;
        b_lb = 0.27;
        E_ub = 10;
        E_lb = 1e-6;
        n_ub = 10;
        n_lb = 1e-2;
        mix_ub=0.5;
        sig_ub = 2500;


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
        % lb_HL(5:5:5*num_sub-5)=mix_lb; % Mixture parameter
        
        A_HL=zeros(1,length(ub_HL));
        A_HL(5:5:5*num_sub-5) = 1;
        b_HL = 1 - lower_mixture_limit;

        ConcT = 1;
        TimeT = 21;


    %% BD initial points
        lb_BD=zeros(11,1);
        ub_BD=lb_BD;
        ub_BD(1)=0.5;
        ub_BD(2)=1;ub_BD(7)=1; % natural birth rate
        % lb_BD(2)=9;
        ub_BD(3)=1;ub_BD(8)=1; % natural death rate
        % lb_BD(3)=9;
        ub_BD(4)=1-1e-6;ub_BD(9)=1-1e-6; % b
        lb_BD(4)=0.27;lb_BD(9)=0.27;
        ub_BD(5)=10;ub_BD(10)=10; % E
        lb_BD(5)=1e-6;lb_BD(10)=1e-6;
        ub_BD(6)=10;ub_BD(11)=10; % n
        lb_BD(6)=1e-2;lb_BD(11)=1e-2;
        lb_BD   = [lb_BD;0];
        ub_BD   = [ub_BD;700];
        A_BD    = [0,1,-1,0,0,0,0,0,0,0,0,0;
             0,0,0,0,0,0,1,-1,0,0,0,0;
             0,-1,1,0,0,0,0,0,0,0,0,0;
             0,0,0,0,0,0,-1,1,0,0,0,0];
        b_BD    = [0.1;0.1;0;0];


        x_init_HL    = [];
        for n=1:num_optim
            x02      = rand(length(ub_HL),1).*(ub_HL-lb_HL) + lb_HL;
        %     x02      = [x02(1);x02(3:12)];
            x_init_HL   = [x_init_HL,x02];
        end
        x_init_BD    = [];
        for n=1:num_optim
            x02      = rand(length(ub_BD),1).*(ub_BD-lb_BD) + lb_BD;
        %     x02      = [x02(1);x02(3:12)];
            x02(3) = max(0,x02(2) - rand*0.1);
            x02(8) = max(x02(7)  - rand*0.1,0);
            x_init_BD   = [x_init_BD,x02];
        end



    %% Generation bound
             lb_GE = lb_BD;
             ub_GE = ub_BD;
             lb_GE(5)  = 0.05;
             ub_GE(5)  = 0.1;
             lb_GE(10) = 0.5;
             ub_GE(10) = 2.5;
             lb_GE(1)  = 0.3;
%              ub_GE(1)  = 0.97;
             lb_GE(4)  = 0.8;
             ub_GE(4)  = 0.9;
             lb_GE(9)  = 0.8;
             ub_GE(9)  = 0.9;
             ub_GE(6)  = 5;
             lb_GE(6)  = 1.5;
             ub_GE(11) = 5;
             lb_GE(11) = 1.5;

             NT      = 13; 
             NC      = 11;
             NR      = 20;


opt_hl_hist  = [];
opt_dyn_hist = [];
opt_sto_hist = [];

ci_hl_p    = [];
ci_hl_GR1  = [];
ci_hl_GR2  = [];
ci_dyn_p   = [];
ci_dyn_GR1 = [];
ci_dyn_GR2 = [];
ci_sto_p   = [];
ci_sto_GR1 = [];
ci_sto_GR2 = [];

Conc = 10^(6)*[0 31.25*10^(-9) 62.5*10^(-9) 125*10^(-9) 250*10^(-9) 375*10^(-9) 500*10^(-9) 1.25*10^(-6) 2.5*10^(-6) 3.75*10^(-6) 5*10^(-6)];
Time1 =   3*(0:NT-1);
% 
theta = rand(12,1).*(ub_GE-lb_GE) + lb_GE;
theta(3) = max(0,theta(2) - rand * 0.1);
theta(8) = max(theta(7) - rand * 0.1,0);
NO       = [100,200,300,400,500];  % Noise level.





Prec_p   = [];
Prec_GR1 = [];
Prec_GR2 = [];
ci_hl_p  = [];
ci_hl_GR1 = [];
ci_hl_GR2 = [];
ci_dyn_p = [];
ci_dyn_GR1 = [];
ci_dyn_GR2 = [];
ci_sto_p = [];
ci_sto_GR1 = [];
ci_sto_GR2 = [];
opt_hl_hist = [];
opt_dyn_hist = [];
opt_sto_hist = [];
% tic
for no  = NO

    theta(12) = no;
    indi_ip  = get_indi(theta,Conc(end));
    

    
    
    %% Performing the bootstrapping 
    B_num = 13;
    B_sample =100;

    Boot_hl  = [];
    Boot_dyn = [];
    Boot_sto = [];
    Boot_hl_GR = [];
    Boot_dyn_GR = [];
    Boot_sto_GR = [];
    B_indi   = [];
    Opt_fval = [];
    [DATA1,Data]  = sto_gen_bd(NR,Conc,Time1,initial,theta,num_sub);



    parfor i = 1:B_sample
        
        Data_i =[];
        indi = [];
        for j = 1:B_num
            ri = randi(20);
            indi = [indi,ri];
            Data_i = [Data_i;DATA1(ri,:,:)];
        end
        B_indi    = [B_indi;indi];
        
        
        func_dyn  = @(x) BD_Obj(x,Data_i,Conc,Time1,B_num,2);
        func_hl   = @(x)vectorized_objective_function_two_noise_levels_k_subpopulations(x,2,Data_i,Conc,Time1, ConcT, TimeT, B_num);
        func_sto  = @(x) BD_obj_func(Data_i, Conc, Time1, B_num, x,2);
        
%             func_hl(theta)
%             func_dyn(theta)
%             func_sto(theta)
%             pause


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
        opt_fval = [opt_fval_hl;opt_fval_dyn;opt_fval_sto];
        Opt_fval = [Opt_fval;opt_fval];
        opt_xx_hl  = params2_hl(:,idx_hl);
        opt_xx_dyn = params2_dyn(:,idx_dyn);
        opt_xx_sto = params2_sto(:,idx_sto);
        opt_xx_hl = GR_sort(opt_xx_hl,Conc(end));
        opt_xx_dyn = GR_sort(opt_xx_dyn,Conc(end));
        opt_xx_sto = GR_sort(opt_xx_sto,Conc(end));
        Boot_hl  = [Boot_hl,opt_xx_hl];
        Boot_dyn = [Boot_dyn,opt_xx_dyn];
        Boot_sto = [Boot_sto,opt_xx_sto];
        Indi_xx_hl  = get_indi(opt_xx_hl,Conc(end));
        Indi_xx_dyn = get_indi(opt_xx_dyn,Conc(end));
        Indi_xx_sto = get_indi(opt_xx_sto,Conc(end));
        Boot_hl_GR  = [Boot_hl_GR,Indi_xx_hl(4:5)'];
        Boot_dyn_GR = [Boot_dyn_GR,Indi_xx_dyn(4:5)'];
        Boot_sto_GR = [Boot_sto_GR,Indi_xx_sto(4:5)'];
    end
    % Compute the confidence interval and record the history
    



    ci_hl_p    = prctile(Boot_hl(5,:),[2.5,97.5]);
    ci_hl_GR1  = prctile(Boot_hl_GR(1,:),[2.5,97.5]);
    ci_hl_GR2  = prctile(Boot_hl_GR(2,:),[2.5,97.5]);
    ci_dyn_p   = prctile(Boot_dyn(1,:),[2.5,97.5]);
    ci_dyn_GR1 = prctile(Boot_dyn_GR(1,:),[2.5,97.5]);
    ci_dyn_GR2 = prctile(Boot_dyn_GR(2,:),[2.5,97.5]);
    ci_sto_p   = prctile(Boot_sto(1,:),[2.5,97.5]);
    ci_sto_GR1 = prctile(Boot_sto_GR(1,:),[2.5,97.5]);
    ci_sto_GR2 = prctile(Boot_sto_GR(2,:),[2.5,97.5]);

    opt_hl_hist  = [opt_hl_hist;Boot_hl];
    opt_dyn_hist = [opt_dyn_hist;Boot_dyn];
    opt_sto_hist = [opt_sto_hist;Boot_sto];
    

%%   Compute the precision
        prec_p_hl    = abs(log(mean(Boot_hl(5,:))/indi_ip(1)));
        prec_GR1_hl  = abs(log(mean(Boot_hl_GR(1,:))/indi_ip(4)));
        prec_GR2_hl  = abs(log(mean(Boot_hl_GR(2,:))/indi_ip(5)));
        prec_p_dyn   = abs(log(mean(Boot_dyn(1,:))/indi_ip(1)));
        prec_GR1_dyn = abs(log(mean(Boot_dyn_GR(1,:))/indi_ip(4)));
        prec_GR2_dyn = abs(log(mean(Boot_dyn_GR(2,:))/indi_ip(5)));
        prec_p_sto   = abs(log(mean(Boot_sto(1,:))/indi_ip(1)));
        prec_GR1_sto = abs(log(mean(Boot_sto_GR(1,:))/indi_ip(4)));
        prec_GR2_sto = abs(log(mean(Boot_sto_GR(2,:))/indi_ip(5)));
    
        prec_p = [prec_p_hl; prec_p_dyn; prec_p_sto];
        prec_GR1 = [prec_GR1_hl; prec_GR1_dyn; prec_GR1_sto];
        prec_GR2 = [prec_GR2_hl; prec_GR2_dyn; prec_GR2_sto];

        Prec_p = [Prec_p, prec_p];
        Prec_GR1 = [Prec_GR1,prec_GR1];
        Prec_GR2 = [Prec_GR2,prec_GR2];
end

hl.p     = ci_hl_p;
hl.GR1   = ci_hl_GR1;
hl.GR2   = ci_hl_GR2;
hl.hist  = opt_hl_hist;
dyn.p    = ci_dyn_p;
dyn.GR1  = ci_dyn_GR1;
dyn.GR2  = ci_dyn_GR2;
dyn.hist = opt_dyn_hist;
sto.p    = ci_sto_p;
sto.GR1  = ci_sto_GR1;
sto.GR2  = ci_sto_GR2;
sto.hist = opt_sto_hist;

% t2 = toc


Info.Conc = Conc;
Info.Time = Time1;
Info.NR   = NR;
Info.B_num = B_num;
Info.num_optim = num_optim;
Info.B_sample = B_sample;
Info.param_GE = [ub_GE,lb_GE];
Info.param_BD = [ub_BD,lb_BD];
Info.param_hl = [ub_HL,lb_HL];
Info.theta    = theta;
Info.B_indi = B_indi;
Info.Opt_fval = Opt_fval;


save('Noise_quant111(var_fixed).mat','DATA1','dyn','hl','sto','Info','Prec_p','Prec_GR1','Prec_GR2');








