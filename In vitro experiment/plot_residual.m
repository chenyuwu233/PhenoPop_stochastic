%% Load data


load('Result\BF_21.mat')






%%  HL

init    = DATA(:,1,1);
Growth  = [exp(opt_xx_hl(1)*Time);exp(opt_xx_hl(6)*Time)];
p       = [opt_xx_hl(5),1-opt_xx_hl(5)];
init_p  = round(init*p);

Pred_hl = init_p*Growth;

%%  DYN

init    = DATA(:,1,1);
ng1     = opt_xx_dyn(2) - opt_xx_dyn(3);
ng2     = opt_xx_dyn(7) - opt_xx_dyn(8);
Growth  = [exp(ng1*Time);exp(ng2*Time)];
p       = [opt_xx_dyn(1),1-opt_xx_dyn(1)];
init_p  = round(init*p);

Pred_dyn = init_p*Growth;


%%  STO

init    = DATA(:,1,1);
ng1     = opt_xx_sto(2) - opt_xx_sto(3);
ng2     = opt_xx_sto(7) - opt_xx_sto(8);
Growth  = [exp(ng1*Time);exp(ng2*Time)];
p       = [opt_xx_sto(1),1-opt_xx_sto(1)];
init_p  = round(init*p);

Pred_sto = init_p*Growth;



%% Histogram plot

% Res  = (squeeze(DATA(:,1,:)) - Pred_sto).^2;

% histogram(Res)


%% Scatter plot

Res_hl  = reshape(squeeze(DATA(:,1,:)) - Pred_hl,1,196);
Res_dyn = reshape(squeeze(DATA(:,1,:)) - Pred_dyn,1,196);
Res_sto = reshape(squeeze(DATA(:,1,:)) - Pred_sto,1,196);
% NR = size(Res_sto,1);
Time_vec = repmat(reshape(repmat(Time,14,1),1,196),1,3);
Time_vec(1:196) = Time_vec(1:196)-0.5;
Time_vec(393:end) = Time_vec(393:end)+0.5;
group   = [ones(1,196),2*ones(1,196),3*ones(1,196)];

s1 = scatter(Time_vec(1:196),Res_hl',136,'filled');
hold on
s2 = scatter(Time_vec(197:392),Res_dyn',136,'filled');
s3 = scatter(Time_vec(393:end),Res_sto',136,'filled');
% gscatter(Time_vec,[Res_hl,Res_dyn,Res_sto],group,'bry','...',[52,52,52])

legend([s1 s2 s3],{'PhenoPop','End-points','Live cell iamge'},'Location','northwest')

ax = gca;
ax.FontWeight = 'bold';
ax.FontSize   = 23;
xlabel('Time')
ylabel('Residual')











