load('Result\CI_BF41_500c.mat')


%%
p_hist = [];
GR1_hist = [];
GR2_hist = [];

%% Compare the confidence interval of p

p = [ci_hl_p(2)-ci_hl_p(1),ci_dyn_p(2)-ci_dyn_p(1),ci_sto_p(2)-ci_sto_p(1)]
p_hist = [p_hist;p];



%% Compare the confidence interval of GR1

GR1 = [ci_hl_GR1(2)-ci_hl_GR1(1),ci_dyn_GR1(2)-ci_dyn_GR1(1),ci_sto_GR1(2)-ci_sto_GR1(1)]
GR1_hist = [GR1_hist;GR1];


%% Compare the confidence interval of p

GR2 = [ci_hl_GR2(2)-ci_hl_GR2(1),ci_dyn_GR2(2)-ci_dyn_GR2(1),ci_sto_GR2(2)-ci_sto_GR2(1)]

GR2_hist = [GR2_hist;GR2];
