function expPowerX = apm_expPowerX(mu,sig,k,pdf)
%This function provides a solution for the expectation value of x^k
%for gaussian: pdf='gauss'
%for gamma: pdf='gamma'
%for lognormal: pdf='lognormal'

switch pdf
    case 'gauss'
        %expPowerX = 1/sqrt(pi) * 2^(k/2 - 1) * sig.^(k-1) ...
        %        .* (sqrt(2) * mu .* gamma(1+k/2) .* hypergeom((1-k)/2,3/2,-mu.^2 ./ (2*sig.^2)) ...
        %        + sig .* gamma((1+k)/2).*hypergeom(-k/2,1/2,- mu.^2 ./ (2*sig.^2)));
        
        %E[x^k]        
        imagFac = (-1i*sqrt(2))^k;        
        expPowerX = imagFac * sig.^k .* kummerU(-k/2,1/2, - mu.^2 ./ (2*sig.^2));
        expPowerX = real(expPowerX);
        %E[|x|^^k]
        %expPowerX = 2^(k/2) * gamma((1+k)/2) / sqrt(pi) * sig.^k .* kummer(-k/2,1/2,- mu.^2 ./ (2*sig.^2));
        %expPowerX = 2^(k/2) * gamma((1+k)/2) / sqrt(pi) * sig.^k .* hypergeom(-k/2,1/2,- mu.^2 ./ (2*sig.^2));
    case 'gamma'
        expPowerX = (sig.^2./mu).^k .* exp(gammaln(k + mu.^2./sig.^2) - gammaln(mu.^2./sig.^2));        
    case 'lognormal'
        %Transform mu and sig to lognormal
        %mu_log = log(mu./sqrt(1 + sig.^2 ./ mu.^2));
        %sig_log = sqrt(log(1 + sig.^2 ./ mu.^2));
        [mu_log,sig_log] = apm_transformMeanStdToLogNormalParameters(mu,sig);
        
        %Enjoy the easy lognormal form
        expPowerX = exp(k * mu_log + k^2 * sig_log.^2 / 2); 
    otherwise
        error(['Expectation value of x^k not available for pdf ''' pdf '''']);

end

