function u = apm_calcCdfDose(d, mu, sigma, cdf_model, is_upper, dmin, dmax)
% Returns CDF(d; par1,par2), where CDF uses the model in cdf_model,
% of parameters par1, par2, set such that the pdf of the dose has the 
% specified mean and standard deviation (mu, sigma).
% If is_upper == 1, returns 1-CDF(d; par1,par2)
if nargin<5
    is_upper = 0;
end
% Support of the shifted beta distribution
if nargin<6
    dmin = 0;
end
if nargin<7
    dmax = 1;
end

switch true
    case (strcmpi(cdf_model,'gauss')==1)
        if is_upper == 0
            u = normcdf(d, mu, sigma);
        else 
            u = normcdf(d, mu, sigma, 'upper');
        end 
   
    case (strcmpi(cdf_model,'loggauss')==1)
        [lnMu, lnVar] = apm_transformMeanCovarianceToLogNormalParameters(mu, sigma^2);
        % Use that the lognormal cdf is the same as the normcdf
        % evaluated at log(d)
        if is_upper == 0
            u = normcdf(log(d), lnMu, lnVar); 
        else 
            u = normcdf(log(d), lnMu, lnVar, 'upper');
        end 
    
    case (strcmpi(cdf_model,'beta')==1)
        [alpha, beta] = apm_transformMeanVarianceToBetaParameters(mu, sigma);
        if is_upper == 0
            u = betacdf(d, alpha, beta);
        else 
            u = betacdf(d, alpha, beta, 'upper');
        end
    
    case (strcmpi(cdf_model,'shiftedbeta')==1)
        [mu_prime,sigma_prime] = apm_transformMeanVarianceShiftedBetaToBetaMeanVariance(mu, sigma, dmin, dmax);
        [alpha,beta] = apm_transformMeanVarianceToBetaParameters(mu_prime, sigma_prime);
        d_prime = (d - dmin) ./ (dmax - dmin);
        u = zeros(1,numel(d_prime));
        if is_upper == 0
            u(d_prime<0) = 0;
            u(d_prime>1) = 1;
            u((d_prime>=0) & (d_prime<=1)) = betacdf(d_prime((d_prime>=0) & (d_prime<=1)), alpha,beta);
        else 
            u(d_prime<0) = 1;
            u(d_prime>1) = 0;
            u((d_prime>=0) & (d_prime<=1)) = betacdf(d_prime((d_prime>=0) & (d_prime<=1)), alpha,beta, 'upper');
        end
        
    
    case (strcmpi(cdf_model,'gamma')==1)
        [k,theta] = apm_transformMeanVarianceToGammaParameters(mu,sigma);
        if is_upper == 0
            u = gamcdf(d, k, theta);
        else 
            u = gamcdf(d, k, theta, 'upper');
        end
    
    case (strcmpi(cdf_model,'gumbel')==1)
        [epsilon,alpha] = apm_transformMeanVarianceToGumbelParameters(mu,sigma);
        if is_upper == 0
            u = evcdf(d, epsilon, alpha);
        else 
            u = evcdf(d, epsilon, alpha, 'upper');
        end
    otherwise
        error(['Invalid copula model']);
end
end
