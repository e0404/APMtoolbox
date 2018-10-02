function pdvh = apm_getPDVHfromProbDVH(expDvh,stdDvh,p,distName,trunc)

if nargin < 5
    trunc = false;
end

if nargin < 4
    distName = 'normal';
end

if numel(expDvh) ~= numel(stdDvh)
    error('Inputs #1 and #2 must have same lengths!');
end

nBins = numel(expDvh);

precCorr = 1e-3;
%To deal with sampled dvhs
stdDvh(stdDvh <= precCorr) = precCorr;
expDvh(expDvh >= (1 - precCorr)) = 1-precCorr;
expDvh(expDvh <= 1e-3) = precCorr;

switch distName
    case 'normal'
        pdvh = arrayfun(@(mu,sigma) getNormalQuantile(mu,sigma,p,trunc),expDvh,stdDvh);
    case 'beta'
        [a,b] = apm_transformMeanVarianceToBetaParameters(expDvh,stdDvh);
        pdvh = arrayfun(@(a,b) getBetaQuantile(a,b,p),a,b);
    otherwise
        error(['Function not implemented for Distributions of type ' distName ]);
end
 
end

function pVal = getNormalQuantile(mu,sigma,p,trunc)

dist = makedist('normal','mu',mu,'sigma',sigma);
if trunc
    dist = truncate(dist,0,1);
end
pVal = dist.icdf(p);
end

function pVal = getBetaQuantile(a,b,p)
dist = makedist('beta','a',a,'b',b);
pVal = dist.icdf(p);
end