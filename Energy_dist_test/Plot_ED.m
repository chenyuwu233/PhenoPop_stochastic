load("ED_100000.mat")

bh = boxplot(ED,{'10','20','50','100','500','1000'});
xlabel('Initial total cell count')
ylabel('Energy distance')
set(gca,'FontSize',25,'FontWeight','bold')
set(bh,'LineWidth',3)