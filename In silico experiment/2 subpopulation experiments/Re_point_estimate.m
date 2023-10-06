parpool('local',20)
warning('off','MATLAB:integral:NonFiniteValue')

for idx = 40:70
    load(strcat('Result\CI',num2str(idx),'(var_fixed).mat'))
    
    %% Initialize
    options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',500,'Display','off','algorithm','sqp');
    num_optim= 100;
    num_sub = 2;
    lower_mixture_limit = 0.10;
    initial = 1000;
    
    ub_HL = Info.param_hl(:,1);
    lb_HL = Info.param_hl(:,2);
    
    ub_BD = Info.param_BD(:,1);
    lb_BD = Info.param_BD(:,2);
    
    Conc  = Info.Conc;
    Time  = Info.Time;
    NR    = Info.NR;
    
    ConcT = 1; % Threshold for concentraltion level
    TimeT = 21; % Threshold for time points
    
    
    x_init_HL    = [];
    for n=1:num_optim
        x02      = rand(length(ub_HL),1).*(ub_HL-lb_HL) + lb_HL;
        x_init_HL   = [x_init_HL,x02];
    end
    x_init_BD    = [];
    for n=1:num_optim
        x02      = rand(length(ub_BD),1).*(ub_BD-lb_BD) + lb_BD;
        x02(3) = max(0,x02(2) - rand*0.1);
        x02(8) = max(x02(7)  - rand*0.1,0);
        x_init_BD   = [x_init_BD,x02];
    end
    
    A_HL=zeros(1,length(ub_HL));
    A_HL(5:5:5*num_sub-5) = 1;
    b_HL = 1 - lower_mixture_limit;
    
    A_BD    = [0,1,-1,0,0,0,0,0,0,0,0,0;
         0,0,0,0,0,0,1,-1,0,0,0,0;
         0,-1,1,0,0,0,0,0,0,0,0,0;
         0,0,0,0,0,0,-1,1,0,0,0,0];
    b_BD    = [0.1;0.1;0;0];
    
    
    %% Optimization
    
    func_dyn  = @(x) BD_Obj(x,DATA1,Conc,Time,NR,2);
    func_hl   = @(x)vectorized_objective_function_two_noise_levels_k_subpopulations(x,2,DATA1,Conc,Time, ConcT, TimeT, NR);
    func_sto  = @(x) BD_obj_func(DATA1, Conc, Time, NR, x,2);
    
    
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
    
    
    %% Dynamic variance model
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
    
    
    %% Stochastic Model
    
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
    
    
    
    [opt_fval_hl, idx_hl]   = min(fval2_hl);
    [opt_fval_dyn, idx_dyn] = min(fval2_dyn);
    [opt_fval_sto, idx_sto] = min(fval2_sto);
    opt_fval = [opt_fval_hl;opt_fval_dyn;opt_fval_sto];
    opt_xx_hl  = params2_hl(:,idx_hl);
    opt_xx_dyn = params2_dyn(:,idx_dyn);
    opt_xx_sto = params2_sto(:,idx_sto);
    opt_xx_hl = GR_sort(opt_xx_hl,Conc(end));
    opt_xx_dyn = GR_sort(opt_xx_dyn,Conc(end));
    opt_xx_sto = GR_sort(opt_xx_sto,Conc(end));
    
    
    save(strcat('Result\CI',num2str(idx),'(point_estimate).mat'))
    clear
end
