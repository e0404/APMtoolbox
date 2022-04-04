function [dmin, dmax] = apm_getShiftedBetaRange(mu, sigma, x, vois)
% Get the support of the shifted beta function (heuristic)

% For data
dmin = max(mu - 3.*sigma, 0); % No negative dose
dmax = mu + 3.*sigma;

if nargin>2 % For simulations with perfect correlation
  dmin = max(mu - 1.25.*sigma, 0); % No negative dose
  dmax = mu + 1.25.*sigma;
	% Set the range to [0, 1+] outside the ill area + a buffer zone
	dmin(find(x<vois(1).xL*0.8)) = 0;
	dmin(find(x>vois(1).xU*0.8)) = 0;
	dmax(find(x<vois(1).xL*0.8)) = max(1, mu(find(x<vois(1).xL*0.8)) + 1.25 .* sigma(find(x<vois(1).xL*0.8)));
	dmax(find(x>vois(1).xU*0.8)) = max(1, mu(find(x>vois(1).xU*0.8)) + 1.25 .* sigma(find(x>vois(1).xU*0.8)));
end
end
