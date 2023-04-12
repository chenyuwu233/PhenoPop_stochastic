
%%
p = [];
GR1 = [];
GR2 = [];
indi_p = zeros(1,3);
indi_GR1 = zeros(1,3);
indi_GR2 = zeros(1,3);
nomal_p = zeros(1,3);
nomal_GR1 = zeros(1,3);
nomal_GR2 = zeros(1,3);

%%
for i = 40:69
    name = strcat('CI',num2str(i),'(var_fixed).mat');
    load(name)
    true_param = get_indi(Info.theta,Info.Conc(end));
    if true_param(1) > hl.p(2) || true_param(1) < hl.p(1)
        indi_p(1) = indi_p(1) + 1;
    elseif true_param(1) > dyn.p(2) || true_param(1) < dyn.p(1)
        indi_p(2) = indi_p(2) + 1;
    elseif true_param(1) > sto.p(2) || true_param(1) < sto.p(1)
        indi_p(3) = indi_p(3) + 1;
    end

    if true_param(4) > hl.GR1(2) || true_param(4) < hl.GR1(1)
        indi_GR1(1) = indi_GR1(1) + 1;
    elseif true_param(4) > dyn.GR1(2) || true_param(4) < dyn.GR1(1)
        indi_GR1(2) = indi_GR1(2) + 1;
    elseif true_param(4) > sto.GR1(2) || true_param(4) < sto.GR1(1)
        indi_GR1(3) = indi_GR1(3) + 1;
    end

    if true_param(5) > hl.GR2(2) || true_param(5) < hl.GR2(1)
        indi_GR2(1) = indi_GR2(1) + 1;
    elseif true_param(5) > dyn.GR2(2) || true_param(5) < dyn.GR2(1)
        indi_GR2(2) = indi_GR2(2) + 1;
    elseif true_param(5) > sto.GR2(2) || true_param(5) < sto.GR2(1)
        indi_GR2(3) = indi_GR2(3) + 1;
    end

    pi = [hl.p(2) - hl.p(1);dyn.p(2) - dyn.p(1);sto.p(2) - sto.p(1)]
    GR1i = [hl.GR1(2) - hl.GR1(1);dyn.GR1(2) - dyn.GR1(1);sto.GR1(2) - sto.GR1(1)]
    GR2i = [hl.GR2(2) - hl.GR2(1);dyn.GR2(2) - dyn.GR2(1);sto.GR2(2) - sto.GR2(1)]
    
    p = [p,pi./sum(pi)]
    GR1 = [GR1,GR1i./sum(GR1i)]
    GR2 = [GR2,GR2i./sum(GR2i)]
end

%% Clean

p(:,end) = []
GR1(:,end) = []
GR2(:,end) = []



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

writetable(p_t,'p(30).xlsx');
writetable(GR1_t,'GR1(30).xlsx');
writetable(GR2_t,'GR2(30).xlsx');



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

%%
for i = 41:70
    if i == 113 || i == 126 ||i == 147 ||i == 149
        continue
    end
    name = append('CI',num2str(i),'(var_fixed).mat');
    load(name)
    Prec_p = [];
    Prec_GR1 = [];
    Prec_GR2 = [];
    for k = 1:1
        theta = Info.theta;
%         theta(12) = No(k);
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
        
        prec_p_hl    = max(log(mean(Boot_hl(5,:))/indi_ip(1)),log(indi_ip(1)/mean(Boot_hl(5,:))));
        prec_GR1_hl  = max(log(mean(Boot_hl_GR(1,:))/indi_ip(4)),log(indi_ip(4)/mean(Boot_hl_GR(1,:))));
        prec_GR2_hl  = max(log(mean(Boot_hl_GR(2,:))/indi_ip(5)),log(indi_ip(5)/mean(Boot_hl_GR(2,:))));
        prec_p_dyn   = max(log(mean(Boot_dyn(1,:))/indi_ip(1)),log(indi_ip(1)/mean(Boot_dyn(1,:))));
        prec_GR1_dyn = max(log(mean(Boot_dyn_GR(1,:))/indi_ip(4)),log(indi_ip(4)/mean(Boot_dyn_GR(1,:))));
        prec_GR2_dyn = max(log(mean(Boot_dyn_GR(2,:))/indi_ip(5)),log(indi_ip(5)/mean(Boot_dyn_GR(2,:))));
        prec_p_sto   = max(log(mean(Boot_sto(1,:))/indi_ip(1)),log(indi_ip(1)/mean(Boot_sto(1,:))));
        prec_GR1_sto = max(log(mean(Boot_sto_GR(1,:))/indi_ip(4)),log(indi_ip(4)/mean(Boot_sto_GR(1,:))));
        prec_GR2_sto = max(log(mean(Boot_sto_GR(2,:))/indi_ip(5)),log(indi_ip(5)/mean(Boot_sto_GR(2,:))));



%         prec_p_hl    = max(abs(log(mean(Boot_hl(5,:))/indi_ip(1))),abs(log(indi_ip(1)/mean(Boot_hl(5,:)))));
%         prec_GR1_hl  = max(abs(log(mean(Boot_hl_GR(1,:))/indi_ip(4))),abs(log(indi_ip(4)/mean(Boot_hl_GR(1,:)))));
%         prec_GR2_hl  = max(abs(log(mean(Boot_hl_GR(2,:))/indi_ip(5))),abs(log(indi_ip(5)/mean(Boot_hl_GR(2,:)))));
%         prec_p_dyn   = max(abs(log(mean(Boot_dyn(1,:))/indi_ip(1))),abs(log(indi_ip(1)/mean(Boot_dyn(1,:)))));
%         prec_GR1_dyn = max(abs(log(mean(Boot_dyn_GR(1,:))/indi_ip(4))),abs(log(indi_ip(4)/mean(Boot_dyn_GR(1,:)))));
%         prec_GR2_dyn = max(abs(log(mean(Boot_dyn_GR(2,:))/indi_ip(5))),abs(log(indi_ip(5)/mean(Boot_dyn_GR(2,:)))));
%         prec_p_sto   = max(abs(log(mean(Boot_sto(1,:))/indi_ip(1))),abs(log(indi_ip(1)/mean(Boot_sto(1,:)))));
%         prec_GR1_sto = max(abs(log(mean(Boot_sto_GR(1,:))/indi_ip(4))),abs(log(indi_ip(4)/mean(Boot_sto_GR(1,:)))));
%         prec_GR2_sto = max(abs(log(mean(Boot_sto_GR(2,:))/indi_ip(5))),abs(log(indi_ip(5)/mean(Boot_sto_GR(2,:)))));

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




%%

h_sto = [];
h_dyn = [];
h_hl  = [];
for i = 1:11
    data_sto_i = sto.hist(i,:);
    data_sto_i = (data_sto_i - mean(data_sto_i))./std(data_sto_i);
    h_sto_i = kstest(data_sto_i);
    h_sto = [h_sto,h_sto_i];

    data_dyn_i = dyn.hist(i,:);
    data_dyn_i = (data_dyn_i - mean(data_dyn_i))./std(data_dyn_i);
    h_dyn_i = kstest(data_dyn_i);
    h_dyn = [h_dyn,h_dyn_i];

    data_hl_i = hl.hist(i,:);
    data_hl_i = (data_hl_i - mean(data_hl_i))./std(data_hl_i);
    h_hl_i = kstest(data_hl_i);
    h_hl = [h_hl,h_hl_i];
end





%% Add the quantile width


for i = 40:70
    name = strcat('CI',num2str(i),'(var_fixed).mat');
    load(name)
    hl_indi = [];
    dyn_indi = [];
    sto_indi = [];
    for j = 1:100
        hl_indi = [hl_indi;get_indi(hl.hist(:,j),Info.Conc(end))];
        dyn_indi = [dyn_indi;get_indi(dyn.hist(:,j),Info.Conc(end))];
        sto_indi = [sto_indi;get_indi(sto.hist(:,j),Info.Conc(end))];
    end
    sto.indi = sto_indi;
    hl.indi = hl_indi;
    dyn.indi = dyn_indi;


    hl.quanti_p = prctile(hl_indi(:,1),[2.5,97.5]);
    hl.quanti_GR1 = prctile(hl_indi(:,4),[2.5,97.5]);
    hl.quanti_GR2 = prctile(hl_indi(:,5),[2.5,97.5]);

    dyn.quanti_p = prctile(dyn_indi(:,1),[2.5,97.5]);
    dyn.quanti_GR1 = prctile(dyn_indi(:,4),[2.5,97.5]);
    dyn.quanti_GR2 = prctile(dyn_indi(:,5),[2.5,97.5]);


    sto.quanti_p = prctile(sto_indi(:,1),[2.5,97.5]);
    sto.quanti_GR1 = prctile(sto_indi(:,4),[2.5,97.5]);
    sto.quanti_GR2 = prctile(sto_indi(:,5),[2.5,97.5]);
    save(name)
    clear

end

%% Record the width

%%
p = [];
GR1 = [];
GR2 = [];
p_norm = [];
GR1_norm = [];
GR2_norm = [];
indi_p = zeros(1,3);
indi_GR1 = zeros(1,3);
indi_GR2 = zeros(1,3);

for i = 40:69
    name = strcat('CI',num2str(i),'(var_fixed).mat');
    load(name)
    true_param = get_indi(Info.theta,Info.Conc(end));
    if true_param(1) > hl.quanti_p(2) || true_param(1) < hl.quanti_p(1)
        indi_p(1) = indi_p(1) + 1;
    elseif true_param(1) > dyn.quanti_p(2) || true_param(1) < dyn.quanti_p(1)
        indi_p(2) = indi_p(2) + 1;
    elseif true_param(1) > sto.quanti_p(2) || true_param(1) < sto.quanti_p(1)
        indi_p(3) = indi_p(3) + 1;
    end

    if true_param(4) > hl.quanti_GR1(2) || true_param(4) < hl.quanti_GR1(1)
        indi_GR1(1) = indi_GR1(1) + 1;
    elseif true_param(4) > dyn.quanti_GR1(2) || true_param(4) < dyn.quanti_GR1(1)
        indi_GR1(2) = indi_GR1(2) + 1;
    elseif true_param(4) > sto.quanti_GR1(2) || true_param(4) < sto.quanti_GR1(1)
        indi_GR1(3) = indi_GR1(3) + 1;
    end

    if true_param(5) > hl.quanti_GR2(2) || true_param(5) < hl.quanti_GR2(1)
        indi_GR2(1) = indi_GR2(1) + 1;
    elseif true_param(5) > dyn.quanti_GR2(2) || true_param(5) < dyn.quanti_GR2(1)
        indi_GR2(2) = indi_GR2(2) + 1;
    elseif true_param(5) > sto.quanti_GR2(2) || true_param(5) < sto.quanti_GR2(1)
        indi_GR2(3) = indi_GR2(3) + 1;
    end

    pi = [hl.quanti_p(2) - hl.quanti_p(1);dyn.quanti_p(2) - dyn.quanti_p(1);sto.quanti_p(2) - sto.quanti_p(1)]
    GR1i = [hl.quanti_GR1(2) - hl.quanti_GR1(1);dyn.quanti_GR1(2) - dyn.quanti_GR1(1);sto.quanti_GR1(2) - sto.quanti_GR1(1)]
    GR2i = [hl.quanti_GR2(2) - hl.quanti_GR2(1);dyn.quanti_GR2(2) - dyn.quanti_GR2(1);sto.quanti_GR2(2) - sto.quanti_GR2(1)]
    
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

writetable(p_t,'p_prctil(30).xlsx');
writetable(GR1_t,'GR1_prctil(30).xlsx');
writetable(GR2_t,'GR2_prctil(30).xlsx');


