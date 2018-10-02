function [expSqDev,stdSqDev] = apm_squaredDevProb(expDose,covDose,dParam)

if isscalar(dParam)
    dParam = ones(size(expDose))*dParam;
end

expDevTerm = sum((expDose-dParam).^2);

varDevTerm = sum(diag(covDose));

expSqDev = expDevTerm + varDevTerm;

if nargout == 2
    stdSqDev = sqrt(2* (trace(covDose^2) + (expDose - dParam)'*covDose*(expDose - dParam)));
end

end

