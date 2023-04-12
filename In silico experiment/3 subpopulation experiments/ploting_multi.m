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
% fontsize(ht,25,"pixels")
axt.FontSize = 25;
axt.Colormap = Color;
title('True initial proportion')
% legend(NAME)
axhl   = nexttile;
hhl    = pie(hl_p);
% fontsize(hhl,25,"pixels")
axhl.FontSize = 25;
axhl.Colormap = Color;
% legend(NAME)
title('Phenopop')
axdyn  = nexttile;
hdyn   = pie(dyn_p);
% fontsize(hdyn,25,"pixels")
axdyn.FontSize = 25;
axdyn.Colormap = Color;
% legend(NAME)
title('End-points')
axsto  = nexttile;
hsto   = pie(sto_p);
% fontsize(hsto,25,"pixels")
axsto.FontSize = 25;
axsto.Colormap = Color;
% legend(NAME)
title('Live cell image')

%%
axGR = nexttile;
axGR.Layout.Tile = 5;
axGR.Layout.TileSpan = [1 4];

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

boxplot(Est_GR(1:3,:)',["Phenopop", "End-points", "Live cell image"],"orientation","horizontal",'Symbol','*r')
boxplot(Est_GR(4:6,:)',["Phenopop", "End-points", "Live cell image"],"orientation","horizontal",'Symbol','*b')
boxplot(Est_GR(7:9,:)',["Phenopop", "End-points", "Live cell image"],"orientation","horizontal",'Symbol','*y')
% boxplot(Est_GR(10:12,:)',["Phenopop model", "End Points model", "Live Cell Image model"],"orientation","horizontal",'Symbol','*g')



title("GR_{50} estimation")
% fontsize(gca,26,"pixels")


% 
% hl_yconf  = [0.8 0.8 1.2 1.2];
% dyn_yconf = [1.8 1.8 2.2 2.2];
% sto_yconf = [2.8 2.8 3.2 3.2];





