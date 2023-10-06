%%
Prec_p_hl = [];
Prec_p_dyn = [];
Prec_GR1_hl = [];
Prec_GR1_dyn = [];
Prec_GR2_hl = [];
Prec_GR2_dyn = [];


%%
Prec_p_hl = [Prec_p_hl;Prec_p(1,1:5)]
Prec_p_dyn = [Prec_p_dyn;Prec_p(2,1:5)]
Prec_GR1_hl = [Prec_GR1_hl;Prec_GR1(1,1:5)]
Prec_GR1_dyn = [Prec_GR1_dyn;Prec_GR1(2,1:5)]
Prec_GR2_hl = [Prec_GR2_hl;Prec_GR2(1,1:5)]
Prec_GR2_dyn = [Prec_GR2_dyn;Prec_GR2(2,1:5)]


%% CI width
p_ci = [];
GR1_ci = [];
GR2_ci = [];


%%
for i = 36:250
    if i == 32 || i == 35 ||i == 38 ||i == 86 ||i==168||i==171||i==173||i==187||i==201||i==205||i==213||i==217||i==238||i==242||i==245
        continue
    end
    name = append('Close_GR',num2str(i),'(var_fixed).mat');
    try
        load(name)
        pi = [hl.p(1,2) - hl.p(1,1);dyn.p(1,2) - dyn.p(1,1);sto.p(1,2) - sto.p(1,1)]; % Note: hl.p: different E_r times 
        GR1i = [hl.GR1(1,2) - hl.GR1(1,1);dyn.GR1(1,2) - dyn.GR1(1,1);sto.GR1(1,2) - sto.GR1(1,1)];
        GR2i = [hl.GR2(1,2) - hl.GR2(1,1);dyn.GR2(1,2) - dyn.GR2(1,1);sto.GR2(1,2) - sto.GR2(1,1)];
        p_ci = [p_ci,pi./sum(pi)];
        GR1_ci = [GR1_ci,GR1i./sum(GR1i)];
        GR2_ci = [GR2_ci,GR2i./sum(GR2i)];
    catch
    end


end


%% 

indi = [ones(1,size(p_ci,2)), 2 * ones(1,size(p_ci,2)), 3* ones(1,size(p_ci,2))];

p_a = [p_ci(1,:),p_ci(2,:),p_ci(3,:)];
p_a = [p_a;indi]';
GR1_a = [GR1_ci(1,:),GR1_ci(2,:),GR1_ci(3,:)];
GR1_a = [GR1_a;indi]';
GR2_a = [GR2_ci(1,:),GR2_ci(2,:),GR2_ci(3,:)];
GR2_a = [GR2_a;indi]';


%%
p_t = array2table(p_a);
GR1_t = array2table(GR1_a);
GR2_t = array2table(GR2_a);

writetable(p_t,'1_p_ci(15).xlsx');
writetable(GR1_t,'1_GR1_ci(15).xlsx');
writetable(GR2_t,'1_GR2_ci(15).xlsx');




%% Record the difference precision
Prec_diff_p_hl = [];
Prec_diff_GR1_hl = [];
Prec_diff_GR2_hl = [];
Prec_diff_p_dyn = [];
Prec_diff_GR1_dyn = [];
Prec_diff_GR2_dyn = [];
Prec_diff_p_sto = [];
Prec_diff_GR1_sto = [];
Prec_diff_GR2_sto = [];
IG       = [0.15,0.3,0.45,0.85,2.0];

%%
for i = 35:250
    if i == 32 ||i == 38 ||i == 86 ||i==168||i==171||i==173||i==187||i==201||i==205||i==213||i==217||i==238||i==242||i==245
        continue
    end
    try
    name = append('Result\Close_GR',num2str(i),'(point_estimate).mat');
        load(name)
        Prec_p = [];
        Prec_GR1 = [];
        Prec_GR2 = [];
        for k = 1:5
%             theta = Info.theta;
            theta(10) = IG(k);
            indi_ip  = get_indi(theta,Conc(end));
            indi_hl  = get_indi(opt_hl_hist(:,k),Conc(end));
            indi_dyn = get_indi(opt_dyn_hist(:,k),Conc(end));
            indi_sto = get_indi(opt_sto_hist(:,k),Conc(end));

%             Boot_hl  = hl.hist(11*k-10:11*k,:);
%             Boot_dyn = dyn.hist(12*k-11:12*k,:);
%             Boot_sto = sto.hist(12*k-11:12*k,:);
%             Boot_hl_GR  = [];
%             Boot_dyn_GR = [];
%             Boot_sto_GR = [];
%             for j = 1:100
%                 indi_hl  = get_indi(Boot_hl(:,j),Info.Conc(end));
%                 indi_dyn = get_indi(Boot_dyn(:,j),Info.Conc(end));
%                 indi_sto = get_indi(Boot_sto(:,j),Info.Conc(end));
%                 Boot_hl_GR  = [Boot_hl_GR,indi_hl(4:5)'];
%                 Boot_dyn_GR = [Boot_dyn_GR,indi_dyn(4:5)'];
%                 Boot_sto_GR = [Boot_sto_GR,indi_sto(4:5)'];
%             end
    %         prec_p_hl    = abs(mean(Boot_hl(5,:)) - indi_ip(1))/indi_ip(1);
    %         prec_GR1_hl  = abs(mean(Boot_hl_GR(1,:)) - indi_ip(4))/indi_ip(4);
    %         prec_GR2_hl  = abs(mean(Boot_hl_GR(2,:)) - indi_ip(5))/indi_ip(5);
    %         prec_p_dyn   = abs(mean(Boot_dyn(1,:)) - indi_ip(1))/indi_ip(1);
    %         prec_GR1_dyn = abs(mean(Boot_dyn_GR(1,:)) - indi_ip(4))/indi_ip(4);
    %         prec_GR2_dyn = abs(mean(Boot_dyn_GR(2,:)) - indi_ip(5))/indi_ip(5);
    %         prec_p_sto   = abs(mean(Boot_sto(1,:)) - indi_ip(1))/indi_ip(1);
    %         prec_GR1_sto = abs(mean(Boot_sto_GR(1,:)) - indi_ip(4))/indi_ip(4);
    %         prec_GR2_sto = abs(mean(Boot_sto_GR(2,:)) - indi_ip(5))/indi_ip(5);
            
    %         prec_p_hl    = max(log(mean(Boot_hl(5,:))/indi_ip(1)),log(indi_ip(1)/mean(Boot_hl(5,:))));
    %         prec_GR1_hl  = max(log(mean(Boot_hl_GR(1,:))/indi_ip(4)),log(indi_ip(4)/mean(Boot_hl_GR(1,:))));
    %         prec_GR2_hl  = max(log(mean(Boot_hl_GR(2,:))/indi_ip(5)),log(indi_ip(5)/mean(Boot_hl_GR(2,:))));
    %         prec_p_dyn   = max(log(mean(Boot_dyn(1,:))/indi_ip(1)),log(indi_ip(1)/mean(Boot_dyn(1,:))));
    %         prec_GR1_dyn = max(log(mean(Boot_dyn_GR(1,:))/indi_ip(4)),log(indi_ip(4)/mean(Boot_dyn_GR(1,:))));
    %         prec_GR2_dyn = max(log(mean(Boot_dyn_GR(2,:))/indi_ip(5)),log(indi_ip(5)/mean(Boot_dyn_GR(2,:))));
    %         prec_p_sto   = max(log(mean(Boot_sto(1,:))/indi_ip(1)),log(indi_ip(1)/mean(Boot_sto(1,:))));
    %         prec_GR1_sto = max(log(mean(Boot_sto_GR(1,:))/indi_ip(4)),log(indi_ip(4)/mean(Boot_sto_GR(1,:))));
    %         prec_GR2_sto = max(log(mean(Boot_sto_GR(2,:))/indi_ip(5)),log(indi_ip(5)/mean(Boot_sto_GR(2,:))));
    
%     
%             prec_p_hl    = abs(log(mean(Boot_hl(5,:))/indi_ip(1)));
%             prec_GR1_hl  = abs(log(mean(Boot_hl_GR(1,:))/indi_ip(4)));
%             prec_GR2_hl  = abs(log(mean(Boot_hl_GR(2,:))/indi_ip(5)));
%             prec_p_dyn   = abs(log(mean(Boot_dyn(1,:))/indi_ip(1)));
%             prec_GR1_dyn = abs(log(mean(Boot_dyn_GR(1,:))/indi_ip(4)));
%             prec_GR2_dyn = abs(log(mean(Boot_dyn_GR(2,:))/indi_ip(5)));
%             prec_p_sto   = abs(log(mean(Boot_sto(1,:))/indi_ip(1)));
%             prec_GR1_sto = abs(log(mean(Boot_sto_GR(1,:))/indi_ip(4)));
%             prec_GR2_sto = abs(log(mean(Boot_sto_GR(2,:))/indi_ip(5)));
    

            prec_p_hl    = abs(log(indi_hl(1)/indi_ip(1)));
            prec_GR1_hl  = abs(log(indi_hl(4)/indi_ip(4)));
            prec_GR2_hl  = abs(log(indi_hl(5)/indi_ip(5)));
            prec_p_dyn   = abs(log(indi_dyn(1)/indi_ip(1)));
            prec_GR1_dyn = abs(log(indi_dyn(4)/indi_ip(4)));
            prec_GR2_dyn = abs(log(indi_dyn(5)/indi_ip(5)));
            prec_p_sto   = abs(log(indi_sto(1)/indi_ip(1)));
            prec_GR1_sto = abs(log(indi_sto(4)/indi_ip(4)));
            prec_GR2_sto = abs(log(indi_sto(5)/indi_ip(5)));
            
    %         prec_p = [prec_p_hl-prec_p_sto; prec_p_dyn-prec_p_sto];
    %         prec_GR1 = [prec_GR1_hl-prec_GR1_sto; prec_GR1_dyn-prec_GR1_sto];
    %         prec_GR2 = [prec_GR2_hl-prec_GR2_sto; prec_GR2_dyn-prec_GR2_sto];
            prec_p = [prec_p_hl; prec_p_dyn; prec_p_sto];
            prec_GR1 = [prec_GR1_hl; prec_GR1_dyn; prec_GR1_sto];
            prec_GR2 = [prec_GR2_hl; prec_GR2_dyn; prec_GR2_sto];
        
            Prec_p = [Prec_p, prec_p];
            Prec_GR1 = [Prec_GR1,prec_GR1];
            Prec_GR2 = [Prec_GR2,prec_GR2];
        end
        Prec_diff_p_hl = [Prec_diff_p_hl;Prec_p(1,:)];
        Prec_diff_GR1_hl = [Prec_diff_GR1_hl;Prec_GR1(1,:)];
        Prec_diff_GR2_hl = [Prec_diff_GR2_hl;Prec_GR2(1,:)];
        Prec_diff_p_dyn = [Prec_diff_p_dyn;Prec_p(2,:)];
        Prec_diff_GR1_dyn = [Prec_diff_GR1_dyn;Prec_GR1(2,:)];
        Prec_diff_GR2_dyn = [Prec_diff_GR2_dyn;Prec_GR2(2,:)];
        Prec_diff_p_sto = [Prec_diff_p_sto;Prec_p(3,:)];
        Prec_diff_GR1_sto = [Prec_diff_GR1_sto;Prec_GR1(3,:)];
        Prec_diff_GR2_sto = [Prec_diff_GR2_sto;Prec_GR2(3,:)];
    catch
    end
end




%% Add the quantile width

for i = 36:250
%     if i == 113 || i == 126 ||i == 147 ||i == 149
%         continue
%     end
    name = append('Close_GR',num2str(i),'(var_fixed).mat');
    try
        load(name)
    
        hl_indi = [];
        dyn_indi = [];
        sto_indi = [];
        for j = 1:100
            hl_indi = [hl_indi;get_indi(hl.hist(1:11,j),Info.Conc(end))];
            dyn_indi = [dyn_indi;get_indi(dyn.hist(1:12,j),Info.Conc(end))];
            sto_indi = [sto_indi;get_indi(sto.hist(1:12,j),Info.Conc(end))];
        end
        sto_close.indi = sto_indi;
        hl_close.indi = hl_indi;
        dyn_close.indi = dyn_indi;
    
    
        hl_close.quanti_p = prctile(hl_close.indi(:,1),[2.5,97.5]);
        hl_close.quanti_GR1 = prctile(hl_close.indi(:,4),[2.5,97.5]);
        hl_close.quanti_GR2 = prctile(hl_close.indi(:,5),[2.5,97.5]);
    
        dyn_close.quanti_p = prctile(dyn_close.indi(:,1),[2.5,97.5]);
        dyn_close.quanti_GR1 = prctile(dyn_close.indi(:,4),[2.5,97.5]);
        dyn_close.quanti_GR2 = prctile(dyn_close.indi(:,5),[2.5,97.5]);
    
    
        sto_close.quanti_p = prctile(sto_close.indi(:,1),[2.5,97.5]);
        sto_close.quanti_GR1 = prctile(sto_close.indi(:,4),[2.5,97.5]);
        sto_close.quanti_GR2 = prctile(sto_close.indi(:,5),[2.5,97.5]);
    
        
        save(name)
        clear
    catch
    end
end




%% CI width
p_norm = [];
GR1_norm = [];
GR2_norm = [];
indi_p = zeros(1,3);
indi_GR1 = zeros(1,3);
indi_GR2 = zeros(1,3);

%% Model comparison 10

for i = 36:250
%     if i == 113 || i == 126 ||i == 147 ||i == 149
%         continue
%     end
    name = append('Close_GR',num2str(i),'(var_fixed).mat');
    try
        load(name)
        true_param = get_indi(Info.theta,Info.Conc(end));
        if true_param(1) > hl_close.quanti_p(2) || true_param(1) < hl_close.quanti_p(1)
            indi_p(1) = indi_p(1) + 1;
        elseif true_param(1) > dyn_close.quanti_p(2) || true_param(1) < dyn_close.quanti_p(1)
            indi_p(2) = indi_p(2) + 1;
        elseif true_param(1) > sto_close.quanti_p(2) || true_param(1) < sto_close.quanti_p(1)
            indi_p(3) = indi_p(3) + 1;
        end
    
        if true_param(4) > hl_close.quanti_GR1(2) || true_param(4) < hl_close.quanti_GR1(1)
            indi_GR1(1) = indi_GR1(1) + 1;
        elseif true_param(4) > dyn_close.quanti_GR1(2) || true_param(4) < dyn_close.quanti_GR1(1)
            indi_GR1(2) = indi_GR1(2) + 1;
        elseif true_param(4) > sto_close.quanti_GR1(2) || true_param(4) < sto_close.quanti_GR1(1)
            indi_GR1(3) = indi_GR1(3) + 1;
        end
    
        if true_param(5) > hl_close.quanti_GR2(2) || true_param(5) < hl_close.quanti_GR2(1)
            indi_GR2(1) = indi_GR2(1) + 1;
        elseif true_param(5) > dyn_close.quanti_GR2(2) || true_param(5) < dyn_close.quanti_GR2(1)
            indi_GR2(2) = indi_GR2(2) + 1;
        elseif true_param(5) > sto_close.quanti_GR2(2) || true_param(5) < sto_close.quanti_GR2(1)
            indi_GR2(3) = indi_GR2(3) + 1;
        end
    
    
        pi = [hl_close.quanti_p(2) - hl_close.quanti_p(1);dyn_close.quanti_p(2) - dyn_close.quanti_p(1);sto_close.quanti_p(2) - sto_close.quanti_p(1)];
        GR1i = [hl_close.quanti_GR1(2) - hl_close.quanti_GR1(1);dyn_close.quanti_GR1(2) - dyn_close.quanti_GR1(1);sto_close.quanti_GR1(2) - sto_close.quanti_GR1(1)];
        GR2i = [hl_close.quanti_GR2(2) - hl_close.quanti_GR2(1);dyn_close.quanti_GR2(2) - dyn_close.quanti_GR2(1);sto_close.quanti_GR2(2) - sto_close.quanti_GR2(1)];
        
    
        p_norm = [p_norm,pi./sum(pi)]
        GR1_norm = [GR1_norm,GR1i./sum(GR1i)]
        GR2_norm = [GR2_norm,GR2i./sum(GR2i)]
    catch
    end
end



