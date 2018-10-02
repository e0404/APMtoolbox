function [lnMu,lnSigma] = apm_transformMeanVarianceToLogNormalParameters(mu, sigma)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
lnSigma = sqrt(log(1 + (sigma.^2 ./ mu.^2)));
lnMu = log(mu) - lnSigma.^2 ./ 2;%log(mu.^2 ./ sqrt(sigma.^2 +mu.^2));

end

