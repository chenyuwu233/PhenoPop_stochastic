function ret = get_indi_mono(theta,d_max,cmd)
    if cmd == 0
        GR50_1 = GR_50(theta(1) - theta(2), theta(3), theta(4), theta(5),d_max);
        ret = [theta(1) - theta(2), GR50_1];
    else
        GR50_1 = GR_50(theta(1) , theta(2), theta(3), theta(4),d_max);
        ret = [theta(1), GR50_1];
    end
end