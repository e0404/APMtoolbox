function [p,b] = apm_transformMeanVarianceToGammaParameters(mu,sigma)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    p = mu.^2 ./ sigma.^2;
    b = mu ./ sigma.^2;
end

