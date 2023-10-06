%% Description
%  This is the function that compute the negative log likelihood of Mixture data
%  by using the CLT approximated.
%  DATA:     with the dimension [Replicant, Dosage, Time]
%  Time:     record all the time point for the data.
%  Dosage:   record all the dosage that we use for differnet replicant
%  theta:    [(p,r,nu,b,E,n)_i,sig] 
%  num_sub:  number of sub population
%  sig:      observation noise ~ N(0,sig)
%  We will return the log-likelihood function for the data set.

function ret = BD_obj_func_mono(DATA, Conc, Time, NR, theta,init)
    NC       = length(Conc);
    ret   = 0;
    for j = 1:NC
        d = Conc(j);
        theta_j = [theta(1),theta(2)-log(theta(3)+(1-theta(3))/(1+(d/theta(4))^(theta(5))))];
        sig2 = theta(6);
        for i = 1:NR     
            X_ij = squeeze(DATA(i,j,:));
            indi = isnan(X_ij);
            X_ij(indi) = [];
            if isempty(X_ij) 
                continue
            end
            Tij  = Time;
            Tij(indi) = [];
            p_ij = CLT_BP_mono(X_ij,theta_j,Tij,init,sig2,0);
            
            if p_ij ~=0
                ret = ret - log(p_ij);
            else
                ret = ret + 1e6;
            end
        end
    end
    if ret == 0
        ret = 1e6;
    end
    
%     ret_theta = Theta;
end