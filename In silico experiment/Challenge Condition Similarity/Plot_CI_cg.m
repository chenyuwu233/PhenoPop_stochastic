%% Load the data and select the experiment
%  We have data index range from 31-40, 81-89, 162-220, 231-250 for
%  similar subpopulation drug reaction.
%  These index represents different random seed that generate the data and
%  true parameters. 
%  Under these experiments, we perform different experiments with
%  different E_r. In particular,
%  31-40: E_r in [0.15,0.3,0.45,0.85,2.0,3.0]
%  81-89: E_r in [0.15,0.3,0.45,0.85,2.0]
%  162-220: E_r in [0.15,0.3,0.45,0.85,2.0]
%  231-250: E_r in [0.15,0.3,0.45,0.85,2.0]
%  Please use 'id' and 'cg' to looking for correspondent experiment. 
%
%  E.G. id = 33, no = 1 means experiment 33 with E_r = 0.15.

id  = 33;
cg  = 1;
d_name = strcat('Result\Close_GR', num2str(id),'(var_fixed).mat');
load(d_name)

%% Retreat the data
IG       = [0.15,0.3,0.45,0.85,2.0,3.0];
Conc = Info.Conc;
hl_p   = [];
hl_GR1 = [];
hl_GR2 = [];
dyn_p   = [];
dyn_GR1 = [];
dyn_GR2 = [];
sto_p   = [];
sto_GR1 = [];
sto_GR2 = [];





for i = 1:100
   hl_indi = get_indi(hl.hist(11*cg-10:11*cg,i),Info.Conc(end));
   hl_p    = [hl_p, hl_indi(1)];
   hl_GR1  = [hl_GR1, hl_indi(4)];
   hl_GR2  = [hl_GR2, hl_indi(5)];
   dyn_indi = get_indi(dyn.hist(12*cg-11:12*cg,i),Info.Conc(end));
   dyn_p    = [dyn_p, dyn_indi(1)];
   dyn_GR1  = [dyn_GR1, dyn_indi(4)];
   dyn_GR2  = [dyn_GR2, dyn_indi(5)];
   sto_indi = get_indi(sto.hist(12*cg-11:12*cg,i),Info.Conc(end));
   sto_p    = [sto_p, sto_indi(1)];
   sto_GR1  = [sto_GR1, sto_indi(4)];
   sto_GR2  = [sto_GR2, sto_indi(5)]; 
end

GR1 = [hl_GR1;dyn_GR1;sto_GR1];
GR2 = [hl_GR2;dyn_GR2;sto_GR2];

% GR1 = [dyn_GR1;sto_GR1];
% GR2 = [dyn_GR2;sto_GR2];
p   = [dyn_p;sto_p];

theta = Info.theta;
theta(10) = IG(cg);
True_indi = get_indi(theta,Info.Conc(end));



%%  Set the color

GR1_color = [255 167 112];
GR2_color = [80  231 235];
GR1_color = GR1_color./255;
GR2_color = GR2_color./255;
Color     = [GR1_color;GR2_color];


%%  Plot proportion (pie)

True_p = [True_indi(1), 1-True_indi(1)];
hl_p   = [mean(hl.p(cg,:)), 1 - mean(hl.p(cg,:))];
dyn_p  = [mean(dyn.p(cg,:)), 1 - mean(dyn.p(cg,:))];
sto_p  = [mean(sto.p(cg,:)), 1 - mean(sto.p(cg,:))];

% label  = {'ST 1: '; 'ST 2: '};

t = tiledlayout('flow');
ax1 = nexttile;
% ax = gca();
h  = pie(True_p);
fontsize(h,26,"pixels")
ax1.Colormap = Color;
% pText = findobj(h,'Type','text');
% pValue = get(pText,'String');
% cText = strcat(label,pValue);
% pText(1).String = cText(1);
% pText(2).String = cText(2);
title("True initial proportion")
ax2 = nexttile;
h  = pie(hl_p);
fontsize(h,26,"pixels")
ax2.Colormap = Color;
% % pText = findobj(h,'Type','text');
% % pValue = get(pText,'String');
% % cText = strcat(label,pValue);
% % pText(1).String = cText(1);
% % pText(2).String = cText(2);
title("PhenoPop")
ax3 = nexttile;
h  = pie(dyn_p);
fontsize(h,26,"pixels")
ax3.Colormap = Color;
% pText = findobj(h,'Type','text');
% pValue = get(pText,'String');
% cText = strcat(label,pValue);
% pText(1).String = cText(1);
% pText(2).String = cText(2);
title("End-points")
ax4 = nexttile;
h  = pie(sto_p);
fontsize(h,26,"pixels")
ax4.Colormap = Color;
% pText = findobj(h,'Type','text');
% pValue = get(pText,'String');
% cText = strcat(label,pValue);
% pText(1).String = cText(1);
% pText(2).String = cText(2);
title("Live cell image")

% gax = geoaxes(t);
% gax.Layout.Tile = 4;
% gax.Layout.TileSpan = [1 3];


%%  Plot Proportion CI

% t = tiledlayout(1,3);
% 
% ax1 = nexttile;
% boxplot(p',["End Points model","Live Cell Image model"]);
% hold on
% yline(True_indi(1),'-.',{'True initial proportion of sub-type 1'});
% xlim([0.5 2.5])
% title("Initial proportion of sub-type 1 estimation")



%%  Plot GR
ax5 = nexttile;
hold on
y_lim = [0 3];
xlim([0.01 5])
ylim(y_lim)
xticks(Conc)



set(gca,"Xscale","log")
xline([True_indi(4) True_indi(5)],'-.')


for i = 1:length(Conc)
    if Conc(i) < True_indi(4) && Conc(i+1)> True_indi(4)
        GR1_int = [Conc(i),Conc(i+1)];
    end
    if Conc(i) < True_indi(5) && Conc(i+1)>True_indi(5)
        GR2_int = [Conc(i),Conc(i+1)];
    end
end

GR1_xconf = [GR1_int(1),GR1_int(2),GR1_int(2),GR1_int(1)];
GR2_xconf = [GR2_int(1),GR2_int(2),GR2_int(2),GR2_int(1)];
GR1_yconf = [0,0,4,4];
GR2_yconf = [0,0,4,4];

fill(GR1_xconf, GR1_yconf, GR1_color,'FaceAlpha',0.3)
fill(GR2_xconf, GR2_yconf, GR2_color,'FaceAlpha',0.3)

hl_yconf  = [0.8 0.8 1.2 1.2];
dyn_yconf = [1.8 1.8 2.2 2.2];
sto_yconf = [2.8 2.8 3.2 3.2];

% dyn_yconf = [0.8 0.8 1.2 1.2];
% sto_yconf = [1.8 1.8 2.2 2.2];

% fill(GR1_xconf,hl_yconf,GR1_color)
% fill(GR1_xconf,dyn_yconf,GR1_color)
% fill(GR1_xconf,sto_yconf,GR1_color)
% fill(GR2_xconf,hl_yconf,GR2_color)
% fill(GR2_xconf,dyn_yconf,GR2_color)
% fill(GR2_xconf,sto_yconf,GR2_color)

boxplot(GR1',["PhenoPop", "End-points", "Live cell image"],"orientation","horizontal",'Symbol','*r')
% boxplot(GR1',["End Points model", "Live Cell Image model"],"orientation","horizontal")

boxplot(GR2',["PhenoPop", "End-points", "Live cell image"],"orientation","horizontal",'Symbol','*b')
% boxplot(GR2',[ "End Points model", "Live Cell Image model"],"orientation","horizontal")
xline(Conc)
title("GR_{50} estimation")
fontsize(gca,26,"pixels")
xticklabels(Conc)

%  Pie chart
ax5.Layout.Tile = 5;
ax5.Layout.TileSpan = [1 4];
% ax5.Layout.TileSpan = [1 3];

%  Boxplot
% ax4.Layout.Tile = 2;
% ax4.Layout.TileSpan = [1 2];


ax1.FontSize = 26;
ax2.FontSize = 26;% Pie chart
ax3.FontSize = 26;% Pie chart
ax4.FontSize = 26;% Pie chart
ax5.FontSize = 26;



