function ucm = apm_createUncertaintyModel(spots,sigmaS,sigmaR,correlationModel)
%APM_CREATEUNCERTAINTYMODEL Summary of this function goes here
%   Detailed explanation goes here

n = numel(spots);

if strcmp(correlationModel,'rand')
    C_S = apm_createCovarianceMatrix(n,sigmaS,'random');
    C_R = apm_createCovarianceMatrix(n,sigmaR,'random');
elseif strcmp(correlationModel,'perfect')
    C_S = sigmaS^2*(ones(n)+1e-10*eye(n)); % perfect correlation
    C_R = sigmaR^2*(ones(n)+1e-10*eye(n));
elseif strcmp(correlationModel,'block')
    C1 = (ones(ceil(n/2))+1e-10*eye(ceil(n/2)));
    C2 = (ones(floor(n/2))+1e-10*eye(floor(n/2)));
    C_S  = [sigmaS*C1 zeros(ceil(n/2));zeros(floor(n/2)) sigmaS*C2];
    C_R  = [sigmaR*C1 zeros(ceil(n/2));zeros(floor(n/2)) sigmaR*C2];
elseif strcmp(correlationModel,'uncorrelated')
    C_S = sigmaS^2*eye(n); % uncorrelated
    C_R = sigmaR^2*eye(n); 
elseif strcmp(correlationModel,'toeplitz')
    C_S = toeplitz(logspace(0,-0.1,n)) .* ((sigmaS*ones(n,1)) * (sigmaS*ones(n,1))');
    C_R = toeplitz(logspace(0,-0.1,n)) .* ((sigmaR*ones(n,1)) * (sigmaR*ones(n,1))');
end


ucm.covSys = C_S;
ucm.covRand = C_R;

end

