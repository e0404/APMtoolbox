function dvh = apm_DVH(d,nBins,dMax) 
    
if nargin < 2
    nBins = 100;
end

if nargin < 3 
    dMax = 1.25*max(d);
end

dvh(1,:) = linspace(0,dMax,nBins);

dvh(2,:) = arrayfun(@(dParam) apm_doseVolume(d,dParam),dvh(1,:));
    
