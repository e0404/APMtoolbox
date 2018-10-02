function pdvh = apm_PDVH(expDose,covDose,p,nBins,dMax)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
if nargin < 4
    nBins = 100;
end

if nargin < 5
    dMax = 1.25*max(expDose);
end

pdvh(1,:) = linspace(0,dMax,nBins);

pdvh(2,:) = arrayfun(@(dParam) f_pdv(expDose,covDose,dParam,p),pdvh(1,:));

end

function pdv = f_pdv(expDose,covDose,dParam,p)
    fVal = normcdf(dParam * ones(numel(expDose),1),expDose,sqrt(diag(covDose)));
    fVal = heaviside(1 - p - fVal);
    pdv = 1/numel(expDose) * sum(fVal);
end

