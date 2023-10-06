%% Load the path and data




%% Select the parameter estimation

% % beta_s vs beta_r
% param1 = 2;
% param2 = 7;
% 
% % b_s vs b_r
% param1 = 4;
% param2 = 9;
% 
% % E_s vs E_r
% param1 = 5;
% param2 = 10;

% m_s vs m_r
param1 = 6;
param2 = 11;


sto1_error = [];
sto2_error = [];
param2_val = [];

for idx = 41:70

    name = strcat('Result\CI',num2str(idx),'(point_estimate).mat');
    load(name)
    
    true1  = Info.theta(param1);
    true2  = Info.theta(param2);
    
    sto1   = sto.hist(param1,:)-true1;
    sto2   = sto.hist(param2,:)-true2;
    
    sto1_error = [sto1_error,sto1];
    sto2_error = [sto2_error,sto2];
    param2_val = [param2_val,true2];
    

end

scatter(sto1_error,sto2_error,'filled')
xlabel('Error of m_s','FontWeight','bold')
ylabel('Error of m_r','FontWeight','bold')
title('Joint confidence region between m_s and m_r')




