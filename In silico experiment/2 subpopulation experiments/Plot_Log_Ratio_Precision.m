load('Result\CI_quant(30).mat')


Prec_p   = [Prec_diff_p_dyn,Prec_diff_p_sto,Prec_diff_p_hl];
Prec_GR1 = [Prec_diff_GR1_dyn,Prec_diff_GR1_sto,Prec_diff_GR1_hl];
Prec_GR2 = [Prec_diff_GR2_dyn,Prec_diff_GR2_sto,Prec_diff_GR2_hl,];

t = tiledlayout(1,3);
ax1 = nexttile;
bh = boxplot(Prec_p,{'End-points','Live cell image','PhenoPop'});
set(bh,'LineWidth',3)
ylim([0,0.5])
title('Initial Proportion')
xlabel('method')
ylabel('Absolute Log Ratio')
ax1.FontSize = 22;
ax1.FontWeight = "bold";
ax2 = nexttile;
bh = boxplot(Prec_GR1,{'End-points','Live cell image','PhenoPop'});
set(bh,'LineWidth',3)
ylim([0,0.5])
title('Sensitive GR_{50}')
ylabel('Absolute Log Ratio')
xlabel('method')
ax2.FontSize = 22;
ax2.FontWeight = "bold";
ax3 = nexttile;
bh = boxplot(Prec_GR2,{'End-points','Live cell image','PhenoPop'});
set(bh,'LineWidth',3)
ylim([0,0.5])
title('Resistant GR_{50}')
ylabel('Absolute Log Ratio')
xlabel('method')
ax3.FontSize = 22;
ax3.FontWeight = "bold";