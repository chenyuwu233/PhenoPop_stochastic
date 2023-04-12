load('Result\CI_width(30).mat')
%% p

boxplot(p_norm',name)
ylim([0,1])
xlabel('method')
ylabel('CI width')
U12 = ranksum(p_norm(1,:),p_norm(2,:));
U23 = ranksum(p_norm(2,:),p_norm(3,:));
U13 = ranksum(p_norm(1,:),p_norm(3,:));
U_stat = [U12,U23,U13];
H = sigstar({[1,2],[2,3],[1,3]},U_stat);
% title('Initial proportion confidence interval width (normalized)')
%% GR1

boxplot(GR1_norm',name)
ylim([0,1])
xlabel('method')
ylabel('CI width')
U12 = ranksum(GR1_norm(1,:),GR1_norm(2,:));
U23 = ranksum(GR1_norm(2,:),GR1_norm(3,:));
U13 = ranksum(GR1_norm(1,:),GR1_norm(3,:));
U_stat = [U12,U23,U13];
H = sigstar({[1,2],[2,3],[1,3]},U_stat);
% title('GR1 confidence interval width (normalized)')

%% GR2

boxplot(GR2_norm',name)
ylim([0,1])
xlabel('method')
ylabel('CI width')
U12 = ranksum(GR2_norm(1,:),GR2_norm(2,:));
U23 = ranksum(GR2_norm(2,:),GR2_norm(3,:));
U13 = ranksum(GR2_norm(1,:),GR2_norm(3,:));
U_stat = [U12,U23,U13];
H = sigstar({[1,2],[2,3],[1,3]},U_stat);
% title('GR2 confidence interval width (normalized)')