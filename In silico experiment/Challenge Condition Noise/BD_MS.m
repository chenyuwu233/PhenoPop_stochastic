% This is a function that calculate the expected rate and variance of population at time
% t and dose concentration d given parameter Params
% Params = (r, d, b, E, n)

function [mean, var] = BD_MS(Params, t, d)
    hill   = log(Params(3) + (1-Params(3))/(1 + (d/Params(4))^(Params(5))));
    d      = Params(2) - hill;
    lambda = Params(1) - d;
    mean   = exp(lambda*t);
    var    = (Params(1) + d)/lambda * (exp(2*lambda *t) - exp(lambda * t));
end