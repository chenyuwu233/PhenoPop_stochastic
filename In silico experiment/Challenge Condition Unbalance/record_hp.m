%%
Prec_p_hl = [];
Prec_p_dyn = [];
Prec_GR1_hl = [];
Prec_GR1_dyn = [];
Prec_GR2_hl = [];
Prec_GR2_dyn = [];


%%
for i = 191:238
    if i == 153 || i == 161 || i == 173 || i == 182 || i == 183 || i == 185 || i == 187 || i == 189
        continue
    end
    name = append('High_pop',num2str(i),'(var_fixed).mat');
    load(name)
    Prec_p_hl = [Prec_p_hl;Prec_p(1,:)];
    Prec_p_dyn = [Prec_p_dyn;Prec_p(2,:)];
    Prec_GR1_hl = [Prec_GR1_hl;Prec_GR1(1,:)];
    Prec_GR1_dyn = [Prec_GR1_dyn;Prec_GR1(2,:)];
    Prec_GR2_hl = [Prec_GR2_hl;Prec_GR2(1,:)];
    Prec_GR2_dyn = [Prec_GR2_dyn;Prec_GR2(2,:)];
end

%%
indi = [0.85 * ones(1,20), 0.9 * ones(1,20), 0.95* ones(1,20),0.99*ones(1,20)]';

p_hl = [];
for i = 1:4
    p_hl = [p_hl;Prec_p_hl(:,i)];
end
p_hl = [p_hl,indi];

p_dyn = [];
for i = 1:4
    p_dyn = [p_dyn;Prec_p_dyn(:,i)];
end
p_dyn = [p_dyn,indi];

GR1_hl = [];
for i = 1:4
    GR1_hl = [GR1_hl;Prec_GR1_hl(:,i)];
end
GR1_hl = [GR1_hl,indi];

GR2_hl = [];
for i = 1:4
    GR2_hl = [GR2_hl;Prec_GR2_hl(:,i)];
end
GR2_hl = [GR2_hl,indi];

GR1_dyn = [];
for i = 1:4
    GR1_dyn = [GR1_dyn;Prec_GR1_dyn(:,i)];
end
GR1_dyn = [GR1_dyn,indi];

GR2_dyn = [];
for i = 1:4
    GR2_dyn = [GR2_dyn;Prec_GR2_dyn(:,i)];
end
GR2_dyn = [GR2_dyn,indi];


%%
p_t = array2table(p_a);
GR1_t = array2table(GR1_a);
GR2_t = array2table(GR2_a);

writetable(p_t,'p(30).xlsx');
writetable(GR1_t,'GR1(30).xlsx');
writetable(GR2_t,'GR2(30).xlsx');



%% CI width
p_ci = [];
GR1_ci = [];
GR2_ci = [];



%% Model comparison

for i = 151:190
    if i == 153 || i == 161 || i == 173 || i == 182 || i == 183 || i == 185 || i == 187 || i == 189
        continue
    end
    name = append('High_pop',num2str(i),'(var_fixed).mat');
    load(name)
    pi = [hl.p(6,2) - hl.p(6,1);dyn.p(6,2) - dyn.p(6,1);sto.p(6,2) - sto.p(6,1)];
    GR1i = [hl.GR1(6,2) - hl.GR1(6,1);dyn.GR1(6,2) - dyn.GR1(6,1);sto.GR1(6,2) - sto.GR1(6,1)];
    GR2i = [hl.GR2(6,2) - hl.GR2(6,1);dyn.GR2(6,2) - dyn.GR2(6,1);sto.GR2(6,2) - sto.GR2(6,1)];
    
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

writetable(p_t,'99_p_ci(32).xlsx');
writetable(GR1_t,'99_GR1_ci(32).xlsx');
writetable(GR2_t,'99_GR2_ci(32).xlsx');




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
IP       = [0.85,0.9,0.95,0.99];

%%
for i =191:238
    if i == 153 || i == 161 || i == 173 || i == 182 || i == 183 || i == 185 || i == 187 || i == 189
        continue
    end
    name = append('High_pop',num2str(i),'(var_fixed).mat');
    load(name)
    Prec_p = [];
    Prec_GR1 = [];
    Prec_GR2 = [];
    for k = 1:4
        theta = Info.theta;
        theta(1) = IP(k);
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
