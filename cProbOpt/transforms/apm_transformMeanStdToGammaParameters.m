function [k,theta] = apm_transformMeanStdToGammaParameters(mu,sigma)
% Return the shape and scale parameters of the gamma distribution
% (these are the input parameters of gampdf) 
    k = mu.^2 ./ sigma.^2; %shape
    theta = sigma.^2 ./ mu; %scale
end

