function [alpha,beta] = apm_transformMeanStdToShiftedBetaParameters(mu,sigma,a,b)
% Return the parameters alpha and beta for a shifted beta-distributed variable
% with support in [a,b], if their location and shape parameters are known.

aux = (a .* b - (a+b) .* mu + mu.^2 + sigma.^2) / (sigma.^2 .* (b-a));

alpha = (a - mu) .* aux;
beta = - (b - mu) .* aux;

if (any(alpha<0) | any(beta<0))
    warning('Moment to shape transformation actually not defined!');
end

end

