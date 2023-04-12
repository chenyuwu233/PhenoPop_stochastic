%% Read the data
load('Result\CI_noise_10(30).mat')






%% p

boxplot(p_norm',name)
ylim([0,1])
xlabel('method')
ylabel('CI width')
H = sigstar({[1,2],[2,3],[1,3]},[0.34,0.022,0.0037]); % The stat is obtained from Rmd: CI_analyze

%% GR1

boxplot(GR1_norm',name)
ylim([0,1])
xlabel('method')
ylabel('CI width')
H = sigstar({[1,2],[2,3],[1,3]},[0.97,0.34,0.31]);


%% GR2

boxplot(GR2_norm',name)
ylim([0,1])
xlabel('method')
ylabel('CI width')
H = sigstar({[1,2],[2,3],[1,3]},[0.88,0.018,0.0073]);
