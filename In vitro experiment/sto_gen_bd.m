%% Stochastic model generator
%  This is a data generator for stochastic model. 
%  Input:
%       NR: number of replicant
%       D : vector of dosage
%       T : vector of time
%       initial: initial cell number
%       theta: coefficient for each sub-group of cell. theta = ((p)_{1:k-1},(r,nu,b,E,n)_{1:k},c)
%       num_sub: number of subgroup
%       model:
%           1. For each time interval, the generation process are
%           independent to each other
%           2. Generate the data for the whole time points simutaneously 

function [ret_1,ret] = sto_gen_bd(NR,Conc,Time,initial,theta,num_sub)
    %% Initialize the parameter
    Param  = zeros(5,num_sub);
%     m_var  = theta(6*num_sub + 1);  % measurement noise
    pi     = theta(1:6:6*num_sub-6);
    pi     = [pi',1 - sum(pi)]; 
    for i = 1:num_sub
       Param(:,i) = theta(num_sub + 5 * (i-1) :num_sub + 5*i-1).';
    end
    [~,ND] = size(Conc);
    [~,NT] = size(Time);
    ret    = zeros(NR,ND,NT);
    for i = 1:NR
        for j = 1:ND
            %% Model 1 Independent generation
%             for k = 1:NT
%                 if k == 1
%                     ret(i,j,k) = initial;
%                 else
%                     if ret(i,j,k-1) == 0
%                         ret(i,j,k) = 0;
%                     else
%                         [ret(i,j,k),pi] = B_samplepath(ret(i,j,k-1),Param,D(j),T(k-1),T(k),pi);
%                     end
%                 end
%             end
%             pi = theta(1:num_sub-1);
%             pi = [pi,1 - sum(pi)]; 
            %% Model 2 Whole time line generation
            temp = Sto_samplepath(initial,Param,Conc(j),Time,pi);
            ret(i,j,:) = temp;
        end
    end
    ret_1 = ret;
    ret_1(:,:,2:end) = max(ret_1(:,:,2:end) + round(normrnd(0,theta(6*num_sub), NR,ND,NT-1)),0 );
%     ret_1 = max(ret + round(normrnd(0,theta(6*num_sub), NR,ND,NT)), 0);  % The measurement noise
end