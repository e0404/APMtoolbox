function [alpha,beta] = apm_transformMeanVarianceToBetaParameters(mu,sigma)
%This function gives alpha and beta parameters for a beta-distributed
%variable, if their shape parameter is known.

%alpha = (mu.^2 - mu.^3 - mu.*(sigma.^2))./ sigma.^2;
%beta = (mu - 1) .* (mu.^2 - mu + sigma.^2) ./ sigma.^2;

alpha = ((1 - mu.^2)./sigma.^2 - 1./mu) .* mu.^2;
beta = alpha.*(1./mu - 1);

if any(sigma.^2 > mu .* (1 - mu))
    warning('Moment to shape transformation actually not defined!');
end

end

