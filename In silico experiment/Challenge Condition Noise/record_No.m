%% Prec
Prec_p_hl = [];
Prec_p_dyn = [];
Prec_GR1_hl = [];
Prec_GR1_dyn = [];
Prec_GR2_hl = [];
Prec_GR2_dyn = [];


%%
Prec_p_hl = [Prec_p_hl;Prec_p(1,:)]
Prec_p_dyn = [Prec_p_dyn;Prec_p(2,:)]
Prec_GR1_hl = [Prec_GR1_hl;Prec_GR1(1,:)]
Prec_GR1_dyn = [Prec_GR1_dyn;Prec_GR1(2,:)]
Prec_GR2_hl = [Prec_GR2_hl;Prec_GR2(1,:)]
Prec_GR2_dyn = [Prec_GR2_dyn;Prec_GR2(2,:)]

%%
indi = [100 * ones(1,36), 200 * ones(1,36), 300* ones(1,36),400*ones(1,36),500 *ones(1,36)]';

p_hl = [];
for i = 1:5
    p_hl = [p_hl;Prec_p_hl(:,i)];
end
p_hl = [p_hl,indi];

p_dyn = [];
for i = 1:5
    p_dyn = [p_dyn;Prec_p_dyn(:,i)];
end
p_dyn = [p_dyn,indi];

GR1_hl = [];
for i = 1:5
    GR1_hl = [GR1_hl;Prec_GR1_hl(:,i)];
end
GR1_hl = [GR1_hl,indi];

GR2_hl = [];
for i = 1:5
    GR2_hl = [GR2_hl;Prec_GR2_hl(:,i)];
end
GR2_hl = [GR2_hl,indi];

GR1_dyn = [];
for i = 1:5
    GR1_dyn = [GR1_dyn;Prec_GR1_dyn(:,i)];
end
GR1_dyn = [GR1_dyn,indi];

GR2_dyn = [];
for i = 1:5
    GR2_dyn = [GR2_dyn;Prec_GR2_dyn(:,i)];
end
GR2_dyn = [GR2_dyn,indi];



%% CI width
p_ci = [];
GR1_ci = [];
GR2_ci = [];
indi_p = zeros(1,3);
indi_GR1 = zeros(1,3);
indi_GR2 = zeros(1,3);


%% Model comparison 10

for i = 91:121
    if i == 113 || i == 126 ||i == 147 ||i == 149
        continue
    end
    name = append('sup_Noise_quant',num2str(i),'(var_fixed).mat');
    load(name)
    true_param = get_indi(Info.theta,Info.Conc(end));
    if true_param(1) > hl.p(1,2) || true_param(1) < hl.p(1,1)
        indi_p(1) = indi_p(1) + 1;
    elseif true_param(1) > dyn.p(1,2) || true_param(1) < dyn.p(1,1)
        indi_p(2) = indi_p(2) + 1;
    elseif true_param(1) > sto.p(1,2) || true_param(1) < sto.p(1,1)
        indi_p(3) = indi_p(3) + 1;
    end

    if true_param(4) > hl.GR1(1,2) || true_param(4) < hl.GR1(1,1)
        indi_GR1(1) = indi_GR1(1) + 1;
    elseif true_param(4) > dyn.GR1(1,2) || true_param(4) < dyn.GR1(1,1)
        indi_GR1(2) = indi_GR1(2) + 1;
    elseif true_param(4) > sto.GR1(1,2) || true_param(4) < sto.GR1(1,1)
        indi_GR1(3) = indi_GR1(3) + 1;
    end

    if true_param(5) > hl.GR2(1,2) || true_param(5) < hl.GR2(1,1)
        indi_GR2(1) = indi_GR2(1) + 1;
    elseif true_param(5) > dyn.GR2(1,2) || true_param(5) < dyn.GR2(1,1)
        indi_GR2(2) = indi_GR2(2) + 1;
    elseif true_param(5) > sto.GR2(1,2) || true_param(5) < sto.GR2(1,1)
        indi_GR2(3) = indi_GR2(3) + 1;
    end


    pi = [hl.p(1,2) - hl.p(1,1);dyn.p(1,2) - dyn.p(1,1);sto.p(1,2) - sto.p(1,1)];
    GR1i = [hl.GR1(1,2) - hl.GR1(1,1);dyn.GR1(1,2) - dyn.GR1(1,1);sto.GR1(1,2) - sto.GR1(1,1)];
    GR2i = [hl.GR2(1,2) - hl.GR2(1,1);dyn.GR2(1,2) - dyn.GR2(1,1);sto.GR2(1,2) - sto.GR2(1,1)];
    
    p_ci = [p_ci,pi./sum(pi)];
    GR1_ci = [GR1_ci,GR1i./sum(GR1i)];
    GR2_ci = [GR2_ci,GR2i./sum(GR2i)];
end


%% Model comparison 50

for i = 91:121
    if i == 113 || i == 126 ||i == 147 ||i == 149
        continue
    end
    name = append('sup_Noise_quant',num2str(i),'(var_fixed).mat');
    load(name)

    true_param = get_indi(Info.theta,Info.Conc(end));
    if true_param(1) > hl.p(end,2) || true_param(1) < hl.p(end,1)
        indi_p(1) = indi_p(1) + 1;
    elseif true_param(1) > dyn.p(end,2) || true_param(1) < dyn.p(end,1)
        indi_p(2) = indi_p(2) + 1;
    elseif true_param(1) > sto.p(end,2) || true_param(1) < sto.p(end,1)
        indi_p(3) = indi_p(3) + 1;
    end

    if true_param(4) > hl.GR1(end,2) || true_param(4) < hl.GR1(end,1)
        indi_GR1(1) = indi_GR1(1) + 1;
    elseif true_param(4) > dyn.GR1(end,2) || true_param(4) < dyn.GR1(end,1)
        indi_GR1(2) = indi_GR1(2) + 1;
    elseif true_param(4) > sto.GR1(end,2) || true_param(4) < sto.GR1(end,1)
        indi_GR1(3) = indi_GR1(3) + 1;
    end

    if true_param(5) > hl.GR2(end,2) || true_param(5) < hl.GR2(end,1)
        indi_GR2(1) = indi_GR2(1) + 1;
    elseif true_param(5) > dyn.GR2(end,2) || true_param(5) < dyn.GR2(end,1)
        indi_GR2(2) = indi_GR2(2) + 1;
    elseif true_param(5) > sto.GR2(end,2) || true_param(5) < sto.GR2(end,1)
        indi_GR2(3) = indi_GR2(3) + 1;
    end

    pi = [hl.p(end,2) - hl.p(end,1);dyn.p(end,2) - dyn.p(end,1);sto.p(end,2) - sto.p(end,1)];
    GR1i = [hl.GR1(end,2) - hl.GR1(end,1);dyn.GR1(end,2) - dyn.GR1(end,1);sto.GR1(end,2) - sto.GR1(end,1)];
    GR2i = [hl.GR2(end,2) - hl.GR2(end,1);dyn.GR2(end,2) - dyn.GR2(end,1);sto.GR2(end,2) - sto.GR2(end,1)];
    
    p_ci = [p_ci,pi./sum(pi)];
    GR1_ci = [GR1_ci,GR1i./sum(GR1i)];
    GR2_ci = [GR2_ci,GR2i./sum(GR2i)];
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

writetable(p_t,'50_p_ci(30).xlsx');
writetable(GR1_t,'50_GR1_ci(30).xlsx');
writetable(GR2_t,'50_GR2_ci(30).xlsx');



%% Noise comparison

p_ci_hl    = [];
p_ci_dyn   = [];
p_ci_sto   = [];
GR1_ci_hl  = [];
GR1_ci_dyn = [];
GR1_ci_sto = [];
GR2_ci_hl  = [];
GR2_ci_dyn = [];
GR2_ci_sto = [];


for i = 111:150
    if i == 113 || i == 126 ||i == 147 ||i == 149
        continue
    end
    name = append('sup_Noise_quant',num2str(i),'(var_fixed).mat');
    load(name)

    pj_hl = [];
    pj_dyn = [];
    pj_sto = [];
    for j = 1:5
        pj_hl = [pj_hl;hl.p(j,2) - hl.p(j,1)];
        pj_dyn = [pj_dyn;dyn.p(j,2) - dyn.p(j,1)];
        pj_sto = [pj_sto;sto.p(j,2) - sto.p(j,1)];
    end
    
    
    p_ci_hl  = [p_ci_hl,pj_hl];
    p_ci_dyn = [p_ci_dyn,pj_dyn];
    p_ci_sto = [p_ci_sto,pj_sto];

    GR1j_hl = [];
    GR1j_dyn = [];
    GR1j_sto = [];
    for j = 1:5
        GR1j_hl = [GR1j_hl;hl.GR1(j,2) - hl.GR1(j,1)];
        GR1j_dyn = [GR1j_dyn;dyn.GR1(j,2) - dyn.GR1(j,1)];
        GR1j_sto = [GR1j_sto;sto.GR1(j,2) - sto.GR1(j,1)];
    end
    
    
    GR1_ci_hl  = [GR1_ci_hl,GR1j_hl];
    GR1_ci_dyn = [GR1_ci_dyn,GR1j_dyn];
    GR1_ci_sto = [GR1_ci_sto,GR1j_sto];


    GR2j_hl = [];
    GR2j_dyn = [];
    GR2j_sto = [];
    for j = 1:5
        GR2j_hl = [GR2j_hl;hl.GR2(j,2) - hl.GR2(j,1)];
        GR2j_dyn = [GR2j_dyn;dyn.GR2(j,2) - dyn.GR2(j,1)];
        GR2j_sto = [GR2j_sto;sto.GR2(j,2) - sto.GR2(j,1)];
    end
    
    
    GR2_ci_hl  = [GR2_ci_hl,GR2j_hl];
    GR2_ci_dyn = [GR2_ci_dyn,GR2j_dyn];
    GR2_ci_sto = [GR2_ci_sto,GR2j_sto];


end



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
No       = [100,200,300,400,500];

%%
for i = 111:142
    if i == 113 || i == 126 ||i == 147 ||i == 149
        continue
    end
    name = append('sup_Noise_quant',num2str(i),'(var_fixed).mat');
    load(name)
    Prec_p = [];
    Prec_GR1 = [];
    Prec_GR2 = [];
    for k = 1:5
        theta = Info.theta;
        theta(12) = No(k);
        indi_ip  = get_indi(theta,Info.Conc(end));
        Boot_hl  = hl.hist(11*k-10:11*k,:);
        Boot_dyn = dyn.hist(12*k-11:12*k,:);
        Boot_sto = sto.hist(12*k-11:12*k,:);
        Boot_hl_GR  = [];
        Boot_dyn_GR = [];
        Boot_sto_GR = [];
        for j = 1:100
            indi_hl  = get_indi(Boot_hl(:,j),Info.Conc(end));
            indi_dyn = get_indi(Boot_dyn(:,j),Info.Conc(end));
            indi_sto = get_indi(Boot_sto(:,j),Info.Conc(end));
            Boot_hl_GR  = [Boot_hl_GR,indi_hl(4:5)'];
            Boot_dyn_GR = [Boot_dyn_GR,indi_dyn(4:5)'];
            Boot_sto_GR = [Boot_sto_GR,indi_sto(4:5)'];
        end
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


        prec_p_hl    = abs(log(mean(Boot_hl(5,:))/indi_ip(1)));
        prec_GR1_hl  = abs(log(mean(Boot_hl_GR(1,:))/indi_ip(4)));
        prec_GR2_hl  = abs(log(mean(Boot_hl_GR(2,:))/indi_ip(5)));
        prec_p_dyn   = abs(log(mean(Boot_dyn(1,:))/indi_ip(1)));
        prec_GR1_dyn = abs(log(mean(Boot_dyn_GR(1,:))/indi_ip(4)));
        prec_GR2_dyn = abs(log(mean(Boot_dyn_GR(2,:))/indi_ip(5)));
        prec_p_sto   = abs(log(mean(Boot_sto(1,:))/indi_ip(1)));
        prec_GR1_sto = abs(log(mean(Boot_sto_GR(1,:))/indi_ip(4)));
        prec_GR2_sto = abs(log(mean(Boot_sto_GR(2,:))/indi_ip(5)));

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
end








%% Add the quantile width

for i = 91:121
    if i == 113 || i == 126 ||i == 147 ||i == 149
        continue
    end
    name = append('sup_Noise_quant',num2str(i),'(var_fixed).mat');
    load(name)

    hl_indi = [];
    dyn_indi = [];
    sto_indi = [];
    for j = 1:100
        hl_indi = [hl_indi;get_indi(hl.hist(1:11,j),Info.Conc(end))];
        dyn_indi = [dyn_indi;get_indi(dyn.hist(1:12,j),Info.Conc(end))];
        sto_indi = [sto_indi;get_indi(sto.hist(1:12,j),Info.Conc(end))];
    end
    sto10.indi = sto_indi;
    hl10.indi = hl_indi;
    dyn10.indi = dyn_indi;

    hl_indi = [];
    dyn_indi = [];
    sto_indi = [];
    for j = 1:100
        hl_indi = [hl_indi;get_indi(hl.hist(end-10:end,j),Info.Conc(end))];
        dyn_indi = [dyn_indi;get_indi(dyn.hist(end-11:end,j),Info.Conc(end))];
        sto_indi = [sto_indi;get_indi(sto.hist(end-11:end,j),Info.Conc(end))];
    end
    sto50.indi = sto_indi;
    hl50.indi = hl_indi;
    dyn50.indi = dyn_indi;

    hl10.quanti_p = prctile(hl10.indi(:,1),[2.5,97.5]);
    hl10.quanti_GR1 = prctile(hl10.indi(:,4),[2.5,97.5]);
    hl10.quanti_GR2 = prctile(hl10.indi(:,5),[2.5,97.5]);

    dyn10.quanti_p = prctile(dyn10.indi(:,1),[2.5,97.5]);
    dyn10.quanti_GR1 = prctile(dyn10.indi(:,4),[2.5,97.5]);
    dyn10.quanti_GR2 = prctile(dyn10.indi(:,5),[2.5,97.5]);


    sto10.quanti_p = prctile(sto10.indi(:,1),[2.5,97.5]);
    sto10.quanti_GR1 = prctile(sto10.indi(:,4),[2.5,97.5]);
    sto10.quanti_GR2 = prctile(sto10.indi(:,5),[2.5,97.5]);


    hl50.quanti_p = prctile(hl50.indi(:,1),[2.5,97.5]);
    hl50.quanti_GR1 = prctile(hl50.indi(:,4),[2.5,97.5]);
    hl50.quanti_GR2 = prctile(hl50.indi(:,5),[2.5,97.5]);

    dyn50.quanti_p = prctile(dyn50.indi(:,1),[2.5,97.5]);
    dyn50.quanti_GR1 = prctile(dyn50.indi(:,4),[2.5,97.5]);
    dyn50.quanti_GR2 = prctile(dyn50.indi(:,5),[2.5,97.5]);


    sto50.quanti_p = prctile(sto50.indi(:,1),[2.5,97.5]);
    sto50.quanti_GR1 = prctile(sto50.indi(:,4),[2.5,97.5]);
    sto50.quanti_GR2 = prctile(sto50.indi(:,5),[2.5,97.5]);
    
    save(name)
    clear
end



%% CI width
p = [];
GR1 = [];
GR2 = [];
p_norm = [];
GR1_norm = [];
GR2_norm = [];
indi_p = zeros(1,3);
indi_GR1 = zeros(1,3);
indi_GR2 = zeros(1,3);

%% Model comparison 10

for i = 91:121
    if i == 113 || i == 126 ||i == 147 ||i == 149
        continue
    end
    name = append('sup_Noise_quant',num2str(i),'(var_fixed).mat');
    load(name)
    true_param = get_indi(Info.theta,Info.Conc(end));
    if true_param(1) > hl10.quanti_p(2) || true_param(1) < hl10.quanti_p(1)
        indi_p(1) = indi_p(1) + 1;
    elseif true_param(1) > dyn10.quanti_p(2) || true_param(1) < dyn10.quanti_p(1)
        indi_p(2) = indi_p(2) + 1;
    elseif true_param(1) > sto10.quanti_p(2) || true_param(1) < sto10.quanti_p(1)
        indi_p(3) = indi_p(3) + 1;
    end

    if true_param(4) > hl10.quanti_GR1(2) || true_param(4) < hl10.quanti_GR1(1)
        indi_GR1(1) = indi_GR1(1) + 1;
    elseif true_param(4) > dyn10.quanti_GR1(2) || true_param(4) < dyn10.quanti_GR1(1)
        indi_GR1(2) = indi_GR1(2) + 1;
    elseif true_param(4) > sto10.quanti_GR1(2) || true_param(4) < sto10.quanti_GR1(1)
        indi_GR1(3) = indi_GR1(3) + 1;
    end

    if true_param(5) > hl10.quanti_GR2(2) || true_param(5) < hl10.quanti_GR2(1)
        indi_GR2(1) = indi_GR2(1) + 1;
    elseif true_param(5) > dyn10.quanti_GR2(2) || true_param(5) < dyn10.quanti_GR2(1)
        indi_GR2(2) = indi_GR2(2) + 1;
    elseif true_param(5) > sto10.quanti_GR2(2) || true_param(5) < sto10.quanti_GR2(1)
        indi_GR2(3) = indi_GR2(3) + 1;
    end


    pi = [hl10.quanti_p(2) - hl10.quanti_p(1);dyn10.quanti_p(2) - dyn10.quanti_p(1);sto10.quanti_p(2) - sto10.quanti_p(1)];
    GR1i = [hl10.quanti_GR1(2) - hl10.quanti_GR1(1);dyn10.quanti_GR1(2) - dyn10.quanti_GR1(1);sto10.quanti_GR1(2) - sto10.quanti_GR1(1)];
    GR2i = [hl10.quanti_GR2(2) - hl10.quanti_GR2(1);dyn10.quanti_GR2(2) - dyn10.quanti_GR2(1);sto10.quanti_GR2(2) - sto10.quanti_GR2(1)];
    
    p = [p,pi]
    GR1 = [GR1,GR1i]
    GR2 = [GR2,GR2i]

    p_norm = [p_norm,pi./sum(pi)]
    GR1_norm = [GR1_norm,GR1i./sum(GR1i)]
    GR2_norm = [GR2_norm,GR2i./sum(GR2i)]
end


%% Model comparison 50

for i = 91:121
    if i == 113 || i == 126 ||i == 147 ||i == 149
        continue
    end
    name = append('sup_Noise_quant',num2str(i),'(var_fixed).mat');
    load(name)
    true_param = get_indi(Info.theta,Info.Conc(end));
    if true_param(1) > hl50.quanti_p(2) || true_param(1) < hl50.quanti_p(1)
        indi_p(1) = indi_p(1) + 1;
    elseif true_param(1) > dyn50.quanti_p(2) || true_param(1) < dyn50.quanti_p(1)
        indi_p(2) = indi_p(2) + 1;
    elseif true_param(1) > sto50.quanti_p(2) || true_param(1) < sto50.quanti_p(1)
        indi_p(3) = indi_p(3) + 1;
    end

    if true_param(4) > hl50.quanti_GR1(2) || true_param(4) < hl50.quanti_GR1(1)
        indi_GR1(1) = indi_GR1(1) + 1;
    elseif true_param(4) > dyn50.quanti_GR1(2) || true_param(4) < dyn50.quanti_GR1(1)
        indi_GR1(2) = indi_GR1(2) + 1;
    elseif true_param(4) > sto50.quanti_GR1(2) || true_param(4) < sto50.quanti_GR1(1)
        indi_GR1(3) = indi_GR1(3) + 1;
    end

    if true_param(5) > hl50.quanti_GR2(2) || true_param(5) < hl50.quanti_GR2(1)
        indi_GR2(1) = indi_GR2(1) + 1;
    elseif true_param(5) > dyn50.quanti_GR2(2) || true_param(5) < dyn50.quanti_GR2(1)
        indi_GR2(2) = indi_GR2(2) + 1;
    elseif true_param(5) > sto50.quanti_GR2(2) || true_param(5) < sto50.quanti_GR2(1)
        indi_GR2(3) = indi_GR2(3) + 1;
    end


    pi = [hl50.quanti_p(2) - hl50.quanti_p(1);dyn50.quanti_p(2) - dyn50.quanti_p(1);sto50.quanti_p(2) - sto50.quanti_p(1)];
    GR1i = [hl50.quanti_GR1(2) - hl50.quanti_GR1(1);dyn50.quanti_GR1(2) - dyn50.quanti_GR1(1);sto50.quanti_GR1(2) - sto50.quanti_GR1(1)];
    GR2i = [hl50.quanti_GR2(2) - hl50.quanti_GR2(1);dyn50.quanti_GR2(2) - dyn50.quanti_GR2(1);sto50.quanti_GR2(2) - sto50.quanti_GR2(1)];
    
    p = [p,pi]
    GR1 = [GR1,GR1i]
    GR2 = [GR2,GR2i]

    p_norm = [p_norm,pi./sum(pi)]
    GR1_norm = [GR1_norm,GR1i./sum(GR1i)]
    GR2_norm = [GR2_norm,GR2i./sum(GR2i)]
end


%%
indi = [ones(1,size(p_norm,2)), 2 * ones(1,size(p_norm,2)), 3* ones(1,size(p_norm,2))];

p_a = [p_norm(1,:),p_norm(2,:),p_norm(3,:)];
p_a = [p_a;indi]';
GR1_a = [GR1_norm(1,:),GR1_norm(2,:),GR1_norm(3,:)];
GR1_a = [GR1_a;indi]';
GR2_a = [GR2_norm(1,:),GR2_norm(2,:),GR2_norm(3,:)];
GR2_a = [GR2_a;indi]';



%%
p_t = array2table(p_a);
GR1_t = array2table(GR1_a);
GR2_t = array2table(GR2_a);

writetable(p_t,'50_p_prctil(30).xlsx');
writetable(GR1_t,'50_GR1_prctil(30).xlsx');
writetable(GR2_t,'50_GR2_prctil(30).xlsx');




