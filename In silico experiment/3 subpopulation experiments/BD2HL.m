function ret = BD2HL(theta,num_sub,low,high)
    ret = zeros(5*num_sub+1,1);
    ret(end-5) = theta(end-5) - theta(end-4);
    ret(end-4) = theta(end-3);
    ret(end-3) = theta(end-2);
    ret(end-2) = theta(end-1);
    ret(end-1) = low;
    ret(end)   = high;
    tp         = 0;
    for i = 1:num_sub-1
        p   = theta(6*i-5);
        tp  = tp+p;
        lam = theta(6*i-4) - theta(6*i-3);
        b   = theta(6*i-2);
        E   = theta(6*i-1);
        n   = theta(6*i);
        ret(5*i-4:5*i) = [lam;b;E;n;p];
    end
end