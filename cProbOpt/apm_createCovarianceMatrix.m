function SIGMA = apm_createCovarianceMatrix(n,sigma,type,rho)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3
    type = 'random';
end

if nargin < 2
    sigma = 1;
end

if size(sigma) == 1
    sigma = ones(n,1)*sigma;
end

switch type
    case 'uncorrelated'
        SIGMA = diag(sigma).^2;
    case 'perfect'
        SIGMA = sigma*sigma';
    case 'random'
        %SIGMA = (sigma*sigma') .* rand(n,n); % generate a random n x n matrix
        rho = 1/sqrt(n) * rand(n,n);
        rho = rho*rho';
        rho(1:(n+1):end) = 1;
        
        SIGMA = sigma*sigma';
        SIGMA = SIGMA.*rho;
    case 'rho'
        if nargin ~= 4
            error('Incorrect number of arguments');
        end
        rho = repmat(rho,n);
        rho(1:(n+1):end) = 1;
        
        SIGMA = sigma*sigma';
        SIGMA = SIGMA.*rho;
    otherwise
        error(['Could not construct covariance matrix of type ''' type '''']);
end

