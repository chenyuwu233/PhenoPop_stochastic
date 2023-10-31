load('Result\CI_multi.mat')
t = tiledlayout('flow');

%%  Set the color

GR1_color = [255 167 112];
GR2_color = [80  231 235];
GR3_color = [244 255 135];
GR4_color = [186 97  255];
GR1_color = GR1_color./255;
GR2_color = GR2_color./255;
GR3_color = GR3_color./255;
GR4_color = GR4_color./255;
Color     = [GR1_color;GR2_color;GR3_color];
%% Plot the proportion
True_p = indi_ip(1:3:end);
hl_p   = mean(Boot_hl_indi(1:3:end,:),2)';
dyn_p  = mean(Boot_dyn_indi(1:3:end,:),2)';
sto_p  = mean(Boot_sto_indi(1:3:end,:),2)';

NAME = [];
for i = 1:num_sub_GE
    name = append('Sub-group ',string(i));
    NAME = [NAME,name];
end
NAME = num2cell(NAME);

axt    = nexttile;

ht     = pie(True_p);
set(ht(2:2:end),'FontSize',25,'FontWeight','bold');
axt.FontSize = 25;
axt.FontWeight = 'bold';
axt.Colormap = Color;

title('True initial proportion')

axhl   = nexttile;

hhl    = pie(hl_p);
set(hhl(2:2:end),'FontSize',25,'FontWeight','bold');
axhl.FontSize = 25;
axhl.FontWeight = 'bold';
axhl.Colormap = Color;

title('PhenoPop')

axdyn  = nexttile;

hdyn   = pie(dyn_p);
set(hdyn(2:2:end),'FontSize',25,'FontWeight','bold');
axdyn.FontSize = 25;
axdyn.FontWeight = 'bold';
axdyn.Colormap = Color;

title('End-points')

axsto  = nexttile;

hsto   = pie(sto_p);
set(hsto(2:2:end),'FontSize',25,'FontWeight','bold');
axsto.FontSize = 25;
axsto.FontWeight = 'bold';
axsto.Colormap = Color;

title('Live cell image')

%%
axGR = nexttile;
axGR.Layout.Tile = 5;
axGR.Layout.TileSpan = [1 4];
axGR.FontSize = 25;
axGR.FontWeight = 'bold';

Est_GR = zeros(3*num_sub_GE,size(Boot_hl_indi,2));
for i = 1:num_sub_GE
    Est_GR(3*i-2,:) = Boot_hl_indi(3*i,:);
    Est_GR(3*i-1,:) = Boot_dyn_indi(3*i,:);
    Est_GR(3*i,:)   = Boot_sto_indi(3*i,:);
end

hold on
y_lim = [0 3];
xlim([0.01 5])
ylim(y_lim)
xticks(Conc)

set(gca,"Xscale","log")
True_GR      = indi_ip(3:3:end);
True_GR_int  = zeros(length(True_GR),2); 
for j = 1:length(True_GR)
    for i = 1:length(Conc)
        if Conc(i) < True_GR(j) && Conc(i+1) > True_GR(j)
            True_GR_int(j,:) = [Conc(i),Conc(i+1)];
        end
    end
end


xline(True_GR,'-.')
xline(Conc)
GR_yconf = [0,0,4,4];
for i = 1:num_sub_GE
    GRi_xconf = [True_GR_int(i,1),True_GR_int(i,2),True_GR_int(i,2),True_GR_int(i,1)];
    fill(GRi_xconf,GR_yconf,Color(i,:),'FaceAlpha',0.3)
end

bh1 = boxplot(Est_GR(1:3,:)',["PhenoPop", "End-points", "Live cell image"],"orientation","horizontal",'Symbol','*r');
bh2 = boxplot(Est_GR(4:6,:)',["PhenoPop", "End-points", "Live cell image"],"orientation","horizontal",'Symbol','*b');
bh3 = boxplot(Est_GR(7:9,:)',["PhenoPop", "End-points", "Live cell image"],"orientation","horizontal",'Symbol','*y');
% boxplot(Est_GR(10:12,:)',["Phenopop model", "End Points model", "Live Cell Image model"],"orientation","horizontal",'Symbol','*g')

set(bh1,'LineWidth',3)
set(bh2,'LineWidth',3)
set(bh3,'LineWidth',3)



title("GR_{50} estimation")
% fontsize(gca,26,"pixels")


% 
% hl_yconf  = [0.8 0.8 1.2 1.2];
% dyn_yconf = [1.8 1.8 2.2 2.2];
% sto_yconf = [2.8 2.8 3.2 3.2];





