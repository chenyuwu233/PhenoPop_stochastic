
load('Result\CI_width(30).mat')
%%  Color

Colormap = [    0.9290    0.6940    0.1250
                0.8500    0.3250    0.0980
                     0    0.4470    0.7410];

%% p(pair)

plot_vecs(p',name,[0,1],'CI width','Initial proportion',Colormap,'signed rank')
xlabel('method')


%% GR1(pair)

plot_vecs(GR1',name,[0,1],'CI width','Sensitive GR_{50}',Colormap,'signed rank')
xlabel('method')

%% GR2(pair)

plot_vecs(GR2',name,[0,1],'CI width','Resistant GR_{50}',Colormap,'signed rank')
xlabel('method')

