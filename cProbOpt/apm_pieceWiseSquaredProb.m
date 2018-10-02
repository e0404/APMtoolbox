function [expPWSq,stdPWSq] = apm_pieceWiseSquaredProb(expDose,covDose,dMin,dMax,p)

if isscalar(dMin)
    dMin = ones(size(expDose))*dMin;
end

if isscalar(dMax)
    dMax = ones(size(expDose))*dMax;
end

if nargin < 5
    p = 1;
end

if isscalar(p)
    p = p*ones(size(expDose));
end

stdDose = sqrt(diag(covDose));
stdDose(stdDose <= 1e-9) = 1e-9;

muTransMax = expDose - dMax;
muTransMin = dMin - expDose;

standardnormcdf = @(x) 0.5*(1 + erf(x/sqrt(2)));
%normpdf = @(x,mu,sig) 1/(sqrt(2*pi)*sig) * exp(-0.5 * ((x - mu) / sig)^2);
zeronormpdf = @(mu,sig) 1/(sqrt(2*pi)*sig) * exp(-0.5 * (mu / sig)^2);


maxSq = sum(p.*arrayfun(@(mu,sig) (mu^2 + sig^2) * standardnormcdf(mu/sig) + mu*sig^2*zeronormpdf(mu,sig),muTransMax,stdDose));
minSq = sum(p.*arrayfun(@(mu,sig) (mu^2 + sig^2) * standardnormcdf(mu/sig) + mu*sig^2*zeronormpdf(mu,sig),muTransMin,stdDose));


expPWSq = maxSq + minSq;

%Not yet derived
stdPWSq = NaN;
end


