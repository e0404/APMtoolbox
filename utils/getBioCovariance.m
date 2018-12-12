function [vMeanAlphaDose,mCovAlphaDose,vMeanSqrtBetaDose,mCovSqrtBetaDose,mCovBioTot]  = getBioCovariance(machine, Spot ,sysBioError,ixTissue,numComp)

vMeanAlphaDose    = zeros(Spot.numOfSpots*numComp,1);
vMeanSqrtBetaDose = zeros(Spot.numOfSpots*numComp,1);
mCovAlphaDose     = zeros(Spot.numOfSpots*numComp);
mCovSqrtBetaDose  = zeros(Spot.numOfSpots*numComp);
mCovBioTot        = zeros(Spot.numOfSpots*numComp*2);

for j = 1:Spot.numOfSpots
   
   ixJ = Spot.ixEnergy(j);
   vMeanAlphaDose(j*numComp-numComp+1:j*numComp)    = machine.data(ixJ).alphaDose(ixTissue).weight;
   vMeanSqrtBetaDose(j*numComp-numComp+1:j*numComp) = machine.data(ixJ).SqrtBetaDose(ixTissue).weight;
    
   for m = 1:Spot.numOfSpots
      
      ixM = Spot.ixEnergy(m); 
      mCovAlphaDose(j*numComp-numComp+1:j*numComp,m*numComp-numComp+1:m*numComp) = ...
         (machine.data(ixJ).alphaDose(ixTissue).weight.*sysBioError) *...
         (machine.data(ixM).alphaDose(ixTissue).weight.*sysBioError)';
      
      mCovSqrtBetaDose(j*numComp-numComp+1:j*numComp,m*numComp-numComp+1:m*numComp) = ...
         (machine.data(ixJ).SqrtBetaDose(ixTissue).weight.*sysBioError) *...
         (machine.data(ixM).SqrtBetaDose(ixTissue).weight.*sysBioError)';

   end
   
end

mCovBioTot(1:Spot.numOfSpots*numComp,1:Spot.numOfSpots*numComp)         = mCovAlphaDose;
mCovBioTot(Spot.numOfSpots*numComp+1:end,Spot.numOfSpots*numComp+1:end) = mCovSqrtBetaDose;

end

