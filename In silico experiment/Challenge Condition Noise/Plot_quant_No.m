load('Result\NO_quant(30_pe).mat')

IP = [100,200,300,400,500];
t = tiledlayout(1,3);
ax1 = nexttile;
hl_p  = plot(IP,mean(Prec_diff_p_hl),'-o');
hold on
dyn_p = plot(IP,mean(Prec_diff_p_dyn),'-o');
sto_p = plot(IP,mean(Prec_diff_p_sto),'-o');
xlabel('Standard deviation of observation noise')
ylabel('Mean Absolute Log Ratio')
yline(0,'-.')
legend([hl_p,dyn_p,sto_p],{'PhenoPop','End-points','Live cell image'},'Location','southeast')
title('Initial proportion estimation')

ax2 = nexttile;
hl_GR1 = plot(IP,mean(Prec_diff_GR1_hl),'-o');
hold on
dyn_GR1 = plot(IP,mean(Prec_diff_GR1_dyn),'-o');
sto_GR1 = plot(IP,mean(Prec_diff_GR1_sto),'-o');
xlabel('Standard deviation of observation noise')
ylabel('Mean Absolute Log Ratio')
yline(0,'-.')
legend([hl_GR1,dyn_GR1,sto_GR1],{'PhenoPop','End-points','Live cell image'},'Location','southeast')
title('Sensitive GR_{50} estimation')

ax3 = nexttile;
hl_GR2 = plot(IP,mean(Prec_diff_GR2_hl),'-o');
hold on
dyn_GR2 = plot(IP,mean(Prec_diff_GR2_dyn),'-o');
sto_GR2 = plot(IP,mean(Prec_diff_GR2_sto),'-o');
xlabel('Standard deviation of observation noise')
ylabel('Mean Absolute Log Ratio')
yline(0,'-.')
legend([hl_GR2,dyn_GR2,sto_GR2],{'PhenoPop','End-points','Live cell image'},'Location','southeast')
title('Resistant GR_{50} estimation')

%% 
