function sampledStdAccRel = apm_calcSampledStdAccRel(n)

Kn = @(n) sqrt((n-1)/2).*exp(gammaln((n-1)/2) - gammaln(n/2));
Vn = @(n) 2*((n-1)/2 - exp(2*gammaln(n/2) - 2*gammaln((n-1)/2)));
sampledStdAccRel = Kn(n).*sqrt(Vn(n))./sqrt(n-1);

end

