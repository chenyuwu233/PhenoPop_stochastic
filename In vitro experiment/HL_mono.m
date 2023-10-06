function [L, N_hi, N_low] = HL_mono(Params,init,DATA,concvec,timevec, ConcT, TimeT, NR)
    %Params are model parameters. DATA is cell counts at each dependent
    %variable condition. size(DATA) = (NR,NC,NT). 
    %concvec is list of concentration considered.
    %timevec is list of times, and NR is number of replicates at each condition.
    %To predict, we multiply DATA(k,j,1) by {p*popfunc(Param1,concvec(j),timevec(i)) + (1-p)*popfunc(Param2,concvec(j),timevec(i))}
    %Predictions are the same for all replicates, but residuals differ between replicates.
    
        %alpha_m = Params(5*m-4);
        %b_m = Params(5*m-3);
        %E_m = Params(5*m-2);
        %n_m = Params(5*m-1);
        %p_m = Params(5*m); %except p_k which is 1 minus sum(p_m)
        sigH=Params(length(Params)-1);
        sigL=Params(length(Params));
        NC=length(concvec);
        NT=length(timevec);
        no_populations=1;
        
        % Initial cell counts to be used in prediction
        time_zero_DATA = init; %size (NR,NC,1)
        rep_time_zero_DATA = repmat(time_zero_DATA,1,1,NT); %size (NR,NC,NT)
        
        % Prediction at time T by exponential growth with rate affected by hill function of concentration, per cell line
        pred_multiplier = vectorized_popfunc(Params(1:4),concvec,timevec);
        %pred_multiplier = p*vectorized_popfunc(Param1,concvec,timevec) ...
        %                    + (1-p)*vectorized_popfunc(Param2,concvec,timevec); % size (1,NC,NT)
        
        repeated_multiplier = repmat(pred_multiplier,NR,1,1); % size (NR,NC,NT)
        
        pred_data = rep_time_zero_DATA .* repeated_multiplier; % size (NR,NC,NT)
        %pred_data = max(0, pred_data); % Should we rectify negative values?
        resid = DATA(:,:,:) - pred_data(:,:,:); % Residuals at time zero do not enter in the calculation (and are zero anyway)
        
        % Count the number of High and Low noise levels
        Hi_Lo_indicators = ones(NR,NC,NT);
        N_hi = sum(Hi_Lo_indicators(:,concvec <= ConcT,timevec(2:NT) >= TimeT), 'all');
        N_low = NR*NC*NT - N_hi;
        
        % Make arrays to divide squared residuals elementwise by
        sig_Hi_Lo_array = sigL*ones(NR,NC,NT);
        sig_Hi_Lo_array(:,concvec <= ConcT,timevec(2:NT) >= TimeT) = sigH;

        resid_terms = resid.^2 ./ (2 .* sig_Hi_Lo_array.^2);
        resid_sum = sum(resid_terms, 'all', 'omitnan'); % size 1
        
        L = resid_sum + N_hi*log(sqrt(2 * pi)*sigH) + N_low*log(sqrt(2 * pi)*sigL);

end
