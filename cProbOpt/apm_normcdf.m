function F = apm_normcdf(x,mu,sigma)
    F = 0.5*(1 + erf((x - mu)./(sigma*sqrt(2))));
end

