% This is an objective function with Branching process, 
% we specific the birth rate and death rate for each sub-group
% Params = ([p_i,r_i, d_i, b_i, E_i, n_i], c)
%
% PMat_i:
% r_i   = birth rate of group i cell
% d_i   = death rate of group i cell
% b_i   = hill parameter beta of group i cell
% E_i   = hill parameter E of group i cell(inflection point)
% n_i   = hill parameter n of group i cell(slop at inflection point)
%
% c     = variance contribution from independent additive measurement noise
% DATA has dimension NR x NC x NT


function ret = BD_Obj(Params,DATA,Conc,Time,NR,n)


% Params = Params.*S; % Scaling the parameter



%% Separate the parameter for each sub-group 


if n == 0
    [NR, NC, NT] = size(DATA);
    num_obs = NR*NT*NC - get_number_discarded_observations(DATA);
    sig=Params(2);
    mean=Params(1);

    L = sum((DATA-mean).^2 ./ (2*sig^2), 'all', 'omitnan') + num_obs*log(sig) ;

    ret = L;

else
    PMat   = zeros(n,5);
    p_vec  = Params(1:6:6*n-6);
    for j = 1:n-1
        PMat(j,:) = Params(j*6-4:j*6);
    end
    PMat(n,:) = Params(6*n-5:6*n-1);
    p_vec     = [p_vec',1 - sum(p_vec)];
    c         = Params(6*n);

    %% Record the dimension

    NC     = length(Conc);
    NT     = length(Time);
    res    = 0;
    sig_log= 0;

    %% Record the residue

    for j = 1:NC
        for k = 1:NR
            for i = 2: NT
                if isnan(DATA(k,j,i))
                    continue
                end
                mu    = [];
                sig2  = [];
                for m = 1:n
                    [mu_n,sig2_n] = BD_MS(PMat(m,:),Time(i),Conc(j));
                    mu     = [mu;mu_n];
                    sig2   = [sig2;sig2_n];
                end
                sig    = DATA(k,j,1) * p_vec * sig2 +  c^2;
                res_ijk= (DATA(k,j,i) - DATA(k,j,1) * p_vec * mu)^2/(2*sig);
                if(isfinite(DATA(k,j,1))&&isfinite(DATA(k,j,i)))
                  res    = res + res_ijk;
                  sig_log= sig_log + log(sqrt(2 *pi * sig));
                end
%                 if sig<=0
%                     sig2
%                     sig
%                     p_vec
%                     pause
%                 end
            end
        end
    end


    ret     = res + sig_log;
end


end