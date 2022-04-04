function [alpha,beta] = apm_transformMeanStdToBetaParameters(mu,sigma)
% Return the parameters alpha and beta for a beta-distributed
% variable, if their shape parameter is known.

alpha = ((mu .* (1 - mu))./ (sigma.^2) - 1) .* mu;
beta = alpha.*(1./mu - 1);

if any(sigma.^2 > mu .* (1 - mu))
    warning('Moment to shape transformation actually not defined!');
end

end

