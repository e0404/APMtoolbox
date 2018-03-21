function [objFunc,gradient] = obFunc(dij,Voxel,w)


dose      = dij.physicalDose * w;
deviation = (dose - Voxel.presDose');

deviation(Voxel.ixNT' == 1 & deviation < 0) = 0;

% calculate objective function value  
objFunc = (Voxel.penalty' .* deviation)' * deviation;


% calculate gradient   
if nargout > 1 
    gradient = 2 * (dij.physicalDose'*(Voxel.penalty'.*deviation));
end



end

