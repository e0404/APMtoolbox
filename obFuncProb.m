function [objFunc,gradient] = obFuncProb(dij,Voxel,mOmega,w)

vOmega    = mOmega * w;
dose      = dij.physicalDoseExp * w;
deviation = (dose - Voxel.presDose');

deviation(Voxel.ixNT' == 1 & deviation < 0) = 0;

% calculate objective function value  
objFunc = (Voxel.penalty' .* deviation)'   * deviation   + (w' * vOmega);


% calculate gradient   
if nargout > 1 
    gradient = 2 * (dij.physicalDoseExp'*(Voxel.penalty'.*deviation)) + vOmega;
end



end

