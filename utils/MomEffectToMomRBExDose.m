function [vMeanRBExD,vStdRBExD ] = MomEffectToMomRBExDose(vMeanEffect, vStdEffect ,vAlpha_x,vBeta_x)

% this functions convertens the first two moments of the biological effect to moments of the RBE weighted dose

vMeanRBExD = zeros(size(vMeanEffect));
vStdRBExD  = zeros(size(vMeanEffect));

gamma = vAlpha_x./(2*vBeta_x);
ix    = gamma > 0 & ~isinf(gamma) ;
x0    = vMeanEffect;
A0    =  sqrt(vBeta_x(ix).^-1 .* x0(ix) + gamma(ix).^2) - gamma(ix);
A1    =  ( 1 * vBeta_x(ix).^-1 )./( 2*sqrt(vBeta_x(ix).^-1 .* x0(ix) + gamma(ix).^2));           % 1 
A2    = -( 1 * vBeta_x(ix).^-2 )./( 4*((vBeta_x(ix).^-1    .* x0(ix) + gamma(ix).^2).^(3/2))) ;  % 2
A4    = -(15 * vBeta_x(ix).^-4 )./(16*((vBeta_x(ix).^-1    .* x0(ix) + gamma(ix).^2).^(7/2))) ;  % 24
    
vMeanRBExD(ix) =   A0 + (1/2)*A2.*vStdEffect(ix).^2 + (1/24)*A4.*vStdEffect(ix).^4;
vStdRBExD(ix)  =   sqrt(A1.^2 .* vStdEffect(ix).^2  +  (0.5 * A2.^2 .* vStdEffect(ix).^4));


end

