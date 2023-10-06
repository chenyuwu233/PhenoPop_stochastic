%% Description:
%  This is a function to obtain the average absolute log ratio accuracy between group of estimation and target value
%
%  Input:
%  - vec: vector of point estimation
%  - aim: target value
%
%  Output:
%  - ret: average absolute log ratio accuracy

function ret = get_log_acc(vec,aim)
    log_rat = log(mean(vec)/aim);
    ret = abs(log_rat);
end