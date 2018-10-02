function dvh = apm_VDH(d,nBins) 
    
if nargin < 2
    nBins = 100;
end
dvh(2,:) = linspace(0,1,nBins);

d = sort(d,'descend');

%dvh(1,:) = arrayfun(@(dParam) apm_%apm_doseVolume(d,dParam),dvh(1,:));
dvh(1,:) = interp1(linspace(0,1,numel(d)),d,dvh(2,:));
    
