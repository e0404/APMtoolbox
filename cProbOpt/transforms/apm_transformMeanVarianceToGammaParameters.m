function [k,theta] = apm_transformMeanVarianceToGammaParameters(mu,sigma)
% Method of the moments for the shape and scale parameters of the gamma distribution
% These are the parameters of gampdf 
    k = mu.^2 ./ sigma.^2; %shape
    theta = sigma.^2 ./ mu; %scale
end

