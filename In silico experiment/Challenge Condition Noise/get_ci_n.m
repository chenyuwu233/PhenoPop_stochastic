function ret = get_ci_n(data)
    pd = fitdist(data,"Normal");
    ci = paramci(pd);
    ret = ci(:,1)';
end