function [mu,sigma] = apm_transformMeanStdToShiftedBetaMeanStd(mu_shifted,sigma_shifted,a,b)
% Return the mean and sigma of the beta distribution, which has support in [0,1], 
% that is equivalent to a shifted beta which has support in [a,b].

mu = (mu_shifted - a) ./ (b - a);
sigma = (sigma_shifted) ./ (b - a);

end

