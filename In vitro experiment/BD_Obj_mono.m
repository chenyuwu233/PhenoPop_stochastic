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


function ret = BD_Obj_mono(Params,DATA,Conc,Time,NR,init)


% Params = Params.*S; % Scaling the parameter



%% Separate the parameter for each sub-group 



    PMat   = Params(1:5);
    c      = Params(6);

    %% Record the dimension

    NC     = length(Conc);
    NT     = length(Time);
    res    = 0;
    sig_log= 0;

    %% Record the residue

    for j = 1:NC
        for k = 1:NR
            for i = 1: NT
                if isnan(DATA(k,j,i))
                    continue
                end
                mu    = [];
                sig2  = [];
                [mu,sig2] = BD_MS(PMat,Time(i),Conc(j));
                sig    = init * sig2 +  c^2;
                res_ijk= (DATA(k,j,i) - init * mu)^2/(2*sig);
                if(isfinite(init)&&isfinite(DATA(k,j,i)))
                  res    = res + res_ijk;
                  sig_log= sig_log + log(sqrt(2 *pi * sig));
                end
                if sig<=0
                    sig2
                    sig
                    p_vec
                    pause
                end
            end
        end
    end


    ret     = res + sig_log;
end
