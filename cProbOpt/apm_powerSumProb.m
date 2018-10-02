function [expPS,stdPS] = apm_powerSumProb(expDose,covDose,k,method)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin < 4
    method = 'int_gauss';
end

stdDose = sqrt(diag(covDose));

expPS = NaN;
stdPS = NaN;

powersum_func = @(d) sum(d.^k);
powersum_grad = @(d) transpose(k * d.^(k - 1));
powersum_hessian = @(d) diag((k-1)*k*d.^(k - 2));


switch method
    case 'int_gauss'
        if k == 1
            expPS = sum(expDose);
            stdPS = sqrt(sum(covDose(:)));
        else
            expPS = apm_expPowerX(expDose,sqrt(diag(covDose)),k,'gauss');
            expPS = sum(real(expPS));
            stdPS = NaN;
        end
    case 'taylor_up'
        [expPS,stdPS] = apm_taylorUncertaintyPropagation(expDose,covDose,powersum_func,powersum_grad);%,powersum_hessian);
        %[expPS,stdPS] = apm_taylorUncertaintyPropagation(expDose,covDose,powersum_func,powersum_grad,powersum_hessian);
        
    case 'int_gamma'
        expPS = apm_expPowerX(expDose,sqrt(diag(covDose)),k,'gamma');
        expPS = sum(expPS);
    case 'int_lognormal'
        %Some functions for logNormal
        %logNormalDist = @(x,mu,sig) exp(-(log(x)-mu).^2 ./ (2*sig^2)) ./ (sqrt(2*pi) * sig * x);
        %logNormalDistTransformed = @(x,mu,sig) logNormalDist(x,log(mu./sqrt(1 + sig.^2./mu.^2)),sqrt(log(1 + sig.^2 ./ mu.^2)));
        %expValueLogInt = @(k,mu,s) exp(k * mu + k^2 * s.^2 / 2);
        
        [mu_log,Sigma_log] = apm_transformMeanCovarianceToLogNormalParameters(expDose,covDose);
        
        expPS = exp(k * mu_log + k^2 * diag(Sigma_log) / 2); 
        expPS = XSum(expPS);
        
        %varPS = 0;
        expTerms = zeros(size(covDose));
        for i=1:numel(expDose)
            for j=1:numel(expDose)
                expTerms(i,j) = k*(mu_log(i) + mu_log(j)) + k^2 / 2 * (Sigma_log(i,i) + Sigma_log(j,j) + 2*Sigma_log(i,j));
                %varPS = varPS + exp(k*(mu_log(i) + mu_log(j)) + k^2 / 2 * (Sigma_log(i,i) + Sigma_log(j,j) + 2*Sigma_log(i,j)));
            end
        end
        
        expTerms = exp(expTerms);
        %expTerms = expTerms - expPS*expPS';
        %expPS = XSum(expPS);
        %maxTerm = max(expTerms(:));
        %expTerms = exp(expTerms - maxTerm);
        %varPS = exp(maxTerm)*XSum(expTerms(:));
        %varPS = XSum(exp(expTerms(:)));
        
        varPS = XSum(expTerms(:));
        
        varPS = varPS - expPS^2;
        stdPS = sqrt(varPS);
        %expPS = expValueLogInt(k,transformMean(expDose,stdDose),transformStd(expDose,stdDose));
        %expPS = apm_expPowerX(expDose,sqrt(diag(covDose)),k,'lognormal');
        
        
    otherwise
        error(['Method ''' method ''' not defined!']);
end

if isnan(stdPS)
    warning(['No analytical form for variance of the power-sum given for method ''' method '''. Using Taylor to approximate variance.']);
        [~,stdPS] = apm_taylorUncertaintyPropagation(expDose,covDose,powersum_func,powersum_grad);%,powersum_hessian);
end
