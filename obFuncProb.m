function [objFunc,gradient] = obFuncProb(dij,Voxel,mOmega,w)

vOmega    = mOmega * w;

if ~isfield(dij,'mAlphaDose')
   d = dij.physicalDoseExp * w;
else
   d = (dij.mAlphaDoseExp * w) + (dij.mSqrtBetaDoseExp * w).^2;
end

deviation = (d - Voxel.presDose');

deviation(Voxel.ixNT' == 1 & deviation < 0) = 0;

% calculate objective function value  
objFunc = (Voxel.penalty' .* deviation)'   * deviation   + (w' * vOmega);

% calculate gradient
delta = 2 * (Voxel.penalty'.*deviation);

if nargout > 1
   if ~isfield(dij,'mAlphaDose')
      gradient  = dij.physicalDoseExp' * delta + 2 * vOmega;
   else
      vBias     = (delta' * dij.mAlphaDoseExp)';
      quadTerm  = dij.mSqrtBetaDoseExp * w;
      mPsi      = (2*(delta.*quadTerm)'* dij.mSqrtBetaDoseExp)';
      gradient  = vBias + mPsi + 2 * vOmega;
   end
end

end

