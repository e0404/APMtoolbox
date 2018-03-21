function PSI_ijlm =  calcSecRangeMom(vLaSi11,vLaSi22,vLaSi12,vLaSi21,Dev_j,Dev_m,w_j,w_m) 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculate the second raw moment. Please note that i
% and l depict voxel indices and j and m pencil beam indices. This function
% allows to calculate the the correlation between spot j and spot m. 
% 
% call
%   PSI_ijlm =  calcSecRangeMom(vLaSi11,vLaSi22,vLaSi12,vLaSi21,Dev_j,Dev_m,w_j,w_m) 
%
% input
%   vLaSi11:        Lambda + Sigma of spot combination j,j 
%   vLaSi22:        Lambda + Sigma of spot combination m,m  
%   vLaSi12:        Lambda + Sigma of spot combination j,m 
%   vLaSi21:        Lambda + Sigma of spot combination m,j 
%   Dev_j:          distance between radiological depth and the Gaussian means of spot j
%   Dev_m:          distance between radiological depth and the Gaussian means of spot m
%   w_j:            weights of Gaussians of spot j
%   w_m:            weights of Gaussians of spot m
%
% output
%   PSI_ijlm:       second raw moment
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018 Hans-Peter Wieser, Niklas Wahl, Philipp Hennig and Mark Bangert
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mW_CovBio = 0;
bioOffset = mW_CovBio + (w_j * w_m');
Det       = vLaSi11*vLaSi22' - (vLaSi12*vLaSi21');
FracDet   = 1./(2*pi*sqrt(Det));
ExpTerm   = FracDet .* exp(-.5./Det.*((Dev_j.^2)*vLaSi22' - bsxfun(@times,(Dev_j*Dev_m'),(vLaSi12+vLaSi21)) + vLaSi11*(Dev_m.^2)'));
PSI_ijlm  = ones(numel(w_j),1)' * (ExpTerm .* bioOffset) * ones(numel(w_m),1);
 
end





