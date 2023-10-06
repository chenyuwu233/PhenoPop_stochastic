load('Result\CI_width_30.mat') % load the confidence interval width result
load('Result\Pe_accuracy_30.mat') % load the point estimation accuracy result

%% Color

Colormap = [    0.9290    0.6940    0.1250
                     0    0.4470    0.7410
                     0    0.4470    0.7410
                0.9290    0.6940    0.1250];

%% p

name = {'PhenoPop(t)','PhenoPop(T)','Live cell image(t)','Live cell image(T)'};
p    = [hl_1_3_p',hl_3_p',sto_1_3_p',sto_3_p'];


bh = boxplot(p,name);
set(gca,'FontSize',23,'FontWeight','bold')
set(bh,'LineWidth',3)
ylim([0,1])
xlabel('method')
ylabel('CI width')
U12 = signrank(hl_1_3_p - hl_3_p);
U34 = signrank(sto_1_3_p - sto_3_p);
U13 = signrank(hl_1_3_p - sto_1_3_p);
U24 = signrank(hl_3_p - sto_3_p);
U_stat = [U12,U34,U13,U24];
H = sigstar({[1,2],[3,4],[1,3],[2,4]},U_stat);
title('Estimation precision of the initial proportion')



% %% GR1
% 
% GR1   = [hl_1_3_GR1',hl_3_GR1',sto_1_3_GR1',sto_3_GR1'];
% 
% bh = boxplot(GR1,name);
% set(gca,'FontSize',23,'FontWeight','bold')
% set(bh,'LineWidth',3)
% ylim([0,1])
% xlabel('method')
% ylabel('CI width')
% U12 = signrank(hl_1_3_GR1 - hl_3_GR1);
% U34 = signrank(sto_1_3_GR1 - sto_3_GR1);
% U13 = signrank(hl_1_3_GR1 - sto_1_3_GR1);
% U24 = signrank(hl_3_GR1 - sto_3_GR1);
% U_stat = [U12,U34,U13,U24];
% H = sigstar({[1,2],[3,4],[1,3],[2,4]},U_stat);
% % title('GR1 confidence interval width (normalized)')
% % 
% % plot_vecs(GR1,name,[0,1],'CI width','Sensitive GR_{50}',Colormap,'signed rank')
% % xlabel('method')
% 
% %% GR2
% 
% bh = boxplot(GR2_norm',name);
% set(gca,'FontSize',23,'FontWeight','bold')
% set(bh,'LineWidth',3)
% ylim([0,1])
% xlabel('method')
% ylabel('CI width')
% U12 = ranksum(GR2_norm(1,:),GR2_norm(2,:));
% U23 = ranksum(GR2_norm(2,:),GR2_norm(3,:));
% U13 = ranksum(GR2_norm(1,:),GR2_norm(3,:));
% U_stat = [U12,U23,U13];
% H = sigstar({[1,2],[2,3],[1,3]},U_stat);
% % title('GR2 confidence interval width (normalized)')