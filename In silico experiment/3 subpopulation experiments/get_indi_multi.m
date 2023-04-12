function ret = get_indi_multi(theta,d_max,num_sub)
    ret = [];
    if length(theta) == 6*num_sub
        p = 0;
        for i = 1:num_sub-1
            theta_i = theta(6*i-5:6*i);
            GR50    = GR_50(theta_i(2) - theta_i(3), theta_i(4),theta_i(5),theta_i(6),d_max);
            ret     = [ret,theta_i(1),theta_i(2) - theta_i(3),GR50];
            p       = p+theta_i(1);
        end
        lam = theta(end-5)-theta(end-4);
        ret = [ret,1-p,lam,GR_50(lam,theta(end-3),theta(end-2),theta(end-1),d_max)];
    else
        p = 0;
        for i = 1:num_sub-1
            theta_i = theta(5*i-4:5*i);
            GR50    = GR_50(theta_i(1),theta_i(2),theta_i(3),theta_i(4),d_max);
            ret     = [ret,theta_i(5),theta_i(1),GR50];
            p       = p+theta_i(5);
        end
        ret = [ret,1-p,theta(end-5),GR_50(theta(end-5),theta(end-4),theta(end-3),theta(end-2),d_max)];
    end
end