name = {'PhenoPop','End-points','Live cell image'};
name_pair = {'PP - EP','PP - LC','EP - LC'};
% load("Result\10_NO_CI_width(30).mat") % Experiment Results with Noise std c = 100
load("Result\50_NO_CI_width(30).mat") % Experiment Results with Noise std c = 500

%% Color

Colormap = [    0.9290    0.6940    0.1250
                0.8500    0.3250    0.0980
                     0    0.4470    0.7410];


%% p

% bh = boxplot(p_norm',name);
% set(gca,'FontSize',23,'FontWeight','bold')
% set(bh,'LineWidth',3)
% ylim([0,1])
% xlabel('method')
% ylabel('CI width')
% U12 = ranksum(p_norm(1,:),p_norm(2,:));
% U23 = ranksum(p_norm(2,:),p_norm(3,:));
% U13 = ranksum(p_norm(1,:),p_norm(3,:));
% U_stat = [U12,U23,U13];
% H = sigstar({[1,2],[2,3],[1,3]},U_stat);
% title('Initial proportion confidence interval width (normalized)')
%% p(pair)

plot_vecs(p',name,[0,1],'CI width','Initial proportion',Colormap,'signed rank')
xlabel('method')

%% GR1

% bh = boxplot(GR1_norm',name);
% set(gca,'FontSize',23,'FontWeight','bold')
% set(bh,'LineWidth',3)
% ylim([0,1])
% xlabel('method')
% ylabel('CI width')
% U12 = ranksum(GR1_norm(1,:),GR1_norm(2,:));
% U23 = ranksum(GR1_norm(2,:),GR1_norm(3,:));
% U13 = ranksum(GR1_norm(1,:),GR1_norm(3,:));
% U_stat = [U12,U23,U13];
% H = sigstar({[1,2],[2,3],[1,3]},U_stat);
% % title('GR1 confidence interval width (normalized)')

%% GR1(pair)

plot_vecs(GR1',name,[0,1],'CI width','Sensitive GR_{50}',Colormap,'signed rank')
xlabel('method')

%% GR2
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


%% GR2(pair)

plot_vecs(GR2',name,[0,1],'CI width','Resistant GR_{50}',Colormap,'signed rank')
xlabel('method')


