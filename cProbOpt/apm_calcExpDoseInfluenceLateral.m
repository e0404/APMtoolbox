function edij = apm_calcExpDoseInfluenceLateral(x,spots,ucm)

nVox = numel(x);
sigma = [spots(:).sigma];
mu = [spots(:).mu];

sigmaAddSqr = sigma.^2 + diag(ucm.covSys)' + diag(ucm.covRand)';

edij = 1./sqrt( 2*pi*ones(nVox,1)*sigmaAddSqr) .* exp( - ( bsxfun(@minus,ones(nVox,1)*mu,x).^2 ) ./ ( 2*ones(nVox,1)*sigmaAddSqr) );

%edij = ones(nVox,1) ./ sqrt( 2*pi*ones(nVox,1)*sigma.^2) .* exp( - ( bsxfun(@minus,ones(nVox,1)*mu,x).^2 ) ./ ( 2*ones(nVox,1)*sigma.^2) );

end

