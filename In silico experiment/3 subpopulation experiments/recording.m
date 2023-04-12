%%   Compute the precision/confidence interval
    prec_hl  = [];
    prec_dyn = [];
    prec_sto = [];
    
    for i = 1:num_sub_GE
        prec_hl_p   = max(log(median(Boot_hl_indi(3*i-2,:))/indi_ip(3*i-2)), log(indi_ip(3*i-2)/median(Boot_hl_indi(3*i-2,:)))); 
        prec_hl_GR  = max(log(median(Boot_hl_indi(3*i,:))/indi_ip(3*i)), log(indi_ip(3*i)/median(Boot_hl_indi(3*i,:))));
        prec_dyn_p  = max(log(median(Boot_dyn_indi(3*i-2,:))/indi_ip(3*i-2)), log(indi_ip(3*i-2)/median(Boot_dyn_indi(3*i-2,:)))); 
        prec_dyn_GR = max(log(median(Boot_dyn_indi(3*i,:))/indi_ip(3*i)), log(indi_ip(3*i)/median(Boot_dyn_indi(3*i,:))));
        prec_sto_p  = max(log(median(Boot_sto_indi(3*i-2,:))/indi_ip(3*i-2)), log(indi_ip(3*i-2)/median(Boot_sto_indi(3*i-2,:)))); 
        prec_sto_GR = max(log(median(Boot_sto_indi(3*i,:))/indi_ip(3*i)), log(indi_ip(3*i)/median(Boot_sto_indi(3*i,:))));
        prec_hl     = [prec_hl,[prec_hl_p;prec_hl_GR]];
        prec_dyn    = [prec_dyn,[prec_dyn_p;prec_dyn_GR]];
        prec_sto    = [prec_sto,[prec_sto_p;prec_sto_GR]];
    end



%%  Cleaning

th = 1e-3;
idx = [];
i = 0;
while i < size(Boot_sto_indi,2)
    i = i+1;
    if Boot_sto_indi(3,i)<th || Boot_hl_indi(3,i)<th ||Boot_dyn_indi(3,i) < th
        idx = [idx,i];
        Boot_sto_indi(:,i) = [];
        Boot_dyn_indi(:,i) = [];
        Boot_hl_indi(:,i) = [];
    end
end




