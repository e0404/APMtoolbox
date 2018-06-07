function dij = apm_calcDoseInfluenceLateral(x,spots)

nVox = numel(x);
sigma = [spots(:).sigma];
mu = [spots(:).mu];

dij = ones(nVox,1) ./ sqrt( 2*pi*ones(nVox,1)*sigma.^2) .* exp( - ( bsxfun(@minus,ones(nVox,1)*mu,x).^2 ) ./ ( 2*ones(nVox,1)*sigma.^2) );

end

