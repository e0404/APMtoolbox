function [lnMu,lnSigma] = apm_transformMeanCovarianceToLogNormalParameters(mu, Sigma)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
lnSigma = real(log(1 + (Sigma ./ (mu*mu'))));
lnSigma(isnan(lnSigma)) = 0;
lnMu = log(mu) - diag(lnSigma) ./ 2;%log(mu.^2 ./ sqrt(sigma.^2 +mu.^2));
lnMu(isinf(lnMu)) = 0;
end

