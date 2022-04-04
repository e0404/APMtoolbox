function [lnMu,lnSigma] = apm_transformMeanStdToLogNormalParameters(mu, sigma)
% Return the parameters of the lognormal distribution

lnSigma = sqrt(log(1 + (sigma.^2 ./ mu.^2)));
lnMu = log(mu) - lnSigma.^2 ./ 2;

end

