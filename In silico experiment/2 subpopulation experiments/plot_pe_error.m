

%% Record the point estimation result

    hl_error  = [];
    dyn_error = [];
    sto_error = [];
    hl_indi_error = [];
    dyn_indi_error = [];
    sto_indi_error = [];
    
    
    
    for idx = 41:70
    
        name = strcat('Result\CI',num2str(idx),'(point_estimate).mat');
        load(name)
        
        theta_PP = [Info.theta(2)-Info.theta(3);Info.theta(4:6);Info.theta(1);
                    Info.theta(7)-Info.theta(8);Info.theta(9:11)];
    
        True_indi = get_indi(Info.theta,max(Info.Conc));
        hl_indi   = get_indi(opt_xx_hl,max(Info.Conc));
        dyn_indi  = get_indi(opt_xx_dyn,max(Info.Conc));
        sto_indi  = get_indi(opt_xx_sto,max(Info.Conc));
    
        hl_error = [hl_error,get_acc(opt_xx_hl(1:end-2),theta_PP)];
        dyn_error = [dyn_error,get_acc(opt_xx_dyn,Info.theta)];
        sto_error = [sto_error,get_acc(opt_xx_sto,Info.theta)];
    
        hl_indi_error = [hl_indi_error;get_log_acc(hl_indi,True_indi)];
        dyn_indi_error = [dyn_indi_error;get_log_acc(dyn_indi,True_indi)];
        sto_indi_error = [sto_indi_error;get_log_acc(sto_indi,True_indi)];
    
    
    
    end
    
    save('Result\PE30.mat','hl_error','dyn_error','sto_error', ...
        'hl_indi_error',"dyn_indi_error","sto_indi_error")


%% Load the result

    load('Result\PE30.mat')


%% Plot the abs error of every parameters

    ax = gca;
    ax.XLim = [-1,12];
    ax.YLim = [-0.3,10];
    ax.FontName = 'Arial';
    xticks([0,1,2,3,4,5,6,7,8,9,10,11])
    xticklabels({'p','\beta_s','\nu_s','b_s','E_s','m_s','\beta_r','\nu_r','b_r','E_r','m_r','c'})
    ax.FontSize = 25;
    ax.FontWeight = 'bold';
    title('Point estimation error')
    
    ax = gca;
    hold on
    
    x = [];
    y = [];
    group = [];
    for i = 1:12 % Error from live cell image
        y = [y,sto_error(i,:)];
        x = [x,linspace(i-1-0.25,i-1+0.25,30)];
        group = [group,ones(1,30)];
    %     s = scatter(ax,x,y,24,'filled','b');
    %     s.MarkerFaceAlpha = 0.5;
    end
    
    for i = 1:12 % Error from end-points
        y = [y,dyn_error(i,:)];
        x = [x,linspace(i-1-0.25,i-1+0.25,30)];
        group = [group,2*ones(1,30)];
    %     s = scatter(ax,x,y,24,'x','r');
    %     s.MarkerFaceAlpha = 0.5;
    end
    
    s = gscatter(x,y,group,'rb','ox',12,'on','Parameters','Relative error');
    % s.MarkerFaceAlpha = 0.5;
    
    yline(0,'LineWidth',3);


%% Plot log ratio error for the important parameters
    Prec_p   = [hl_indi_error(:,1),dyn_indi_error(:,1),sto_indi_error(:,1)];
    Prec_GR1 = [hl_indi_error(:,4),dyn_indi_error(:,4),sto_indi_error(:,4)];
    Prec_GR2 = [hl_indi_error(:,5),dyn_indi_error(:,5),sto_indi_error(:,5)];
    
    name     = {'PhenoPop','End-points','Live cell image'};
    
    
    %%  Color
    
    Colormap = [    0.9290    0.6940    0.1250
                    0.8500    0.3250    0.0980
                         0    0.4470    0.7410];
    
    
    %% Initial Proportion
    plot_vecs(Prec_p,name,[0,0.5],'Absolute log ratio','Initial Proportion',Colormap,'signed rank')
    xlabel('method')
    
    %% Sensitive
    plot_vecs(Prec_GR1,name,[0,0.5],'Absolute log ratio','Sensitive GR_{50}',Colormap,'signed rank')
    xlabel('method')
    
    %% Resistant
    plot_vecs(Prec_GR2,name,[0,0.5],'Absolute log ratio','Resistant GR_{50}',Colormap,'signed rank')
    xlabel('method')

