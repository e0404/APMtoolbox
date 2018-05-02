function pwSq = apm_pieceWiseSquared(dose,dMin,dMax)

if isscalar(dMin)
    dMin = ones(size(dose))*dMin;
end

if isscalar(dMax)
    dMax = ones(size(dose))*dMax;
end

%Transform to center
dMaxDev = dose - dMax;
dMinDev = dMin - dose;

%Apply positivity operator
dMaxDev(dMaxDev < 0) = 0;
dMinDev(dMinDev < 0) = 0;

pwSq = sum(arrayfun(@(dMinDev,dMaxDev) dMinDev^2 + dMaxDev^2,dMinDev,dMaxDev));

end

