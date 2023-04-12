% This is a function for sorting the parameter according to the GR-50 from
% small to large 


function ret = GR_sort_multi(par,d_max,num_sub)
    ret = [];
    if length(par) == 6*num_sub
        indi = get_indi_multi(par,d_max,num_sub);
        GR = indi(3:3:end);
        GR_sort = sort(GR);
    elseif length(par) == 5*num_sub+1
        indi = get_indi_multi(par,d_max,num_sub);
        GR = indi(3:3:end);
        GR_sort = sort(GR);
    end
    for i = 1:length(GR_sort)
        idx = find(GR==GR_sort(i));
        ret = [ret,indi(idx*3-2:idx*3)];
    end
end