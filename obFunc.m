function [objFunc,gradient] = obFunc(dij,Voxel,w)

if ~isfield(dij,'mAlphaDose')
   d = dij.physicalDose * w;
else
   d = (dij.mAlphaDose * w) + (dij.mSqrtBetaDose * w).^2;
end

deviation = (d - Voxel.presDose');

deviation(Voxel.ixNT' == 1 & deviation < 0) = 0;

% calculate objective function value
objFunc = (Voxel.penalty' .* deviation)' * deviation;


% calculate gradient
delta = 2 * (Voxel.penalty'.*deviation);

if nargout > 1
   if ~isfield(dij,'mAlphaDose')
      gradient  = dij.physicalDose' * delta;
   else
      vBias     = (delta' * dij.mAlphaDose)';
      quadTerm  = dij.mSqrtBetaDose * w;
      mPsi      = (2*(delta.*quadTerm)'* dij.mSqrtBetaDose)';
      gradient  = vBias + mPsi ; 
   end
end



end

