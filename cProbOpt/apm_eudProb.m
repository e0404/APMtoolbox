function [expEUD,stdEUD] = apm_eudProb(expDose,covDose,k,method)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

expEUD = NaN;
stdEUD = NaN;

fullVoxNum = numel(expDose);


% remove the zero variance voxels
%varDose = diag(covDose);
%ix = find(varDose);
%covDose = covDose(ix,ix);
%expDose = expDose(ix);


% make variance non-zero
% varDose = diag(covDose);
% addDiag = zeros(size(varDose));
% addDiag(varDose == 0) = 1e-9;
% covDose = covDose + diag(addDiag);

%reducedVoxNum = numel(expDose);

%nFac = 1/max(expDose(:));
%expDose = expDose * nFac;
%covDose = covDose * nFac^2;

if nargin < 4
    method = 'taylor_direct';
end

        eud_func = @(d,k) (sum(d.^k)/fullVoxNum)^(1/k);
        eud_grad = @(d,k) transpose(eud_func(d,k)/sum(d.^k) * d.^(k-1));
        %eud_hessian = @(d,k) eud_func(d,k)/sum(d.^k)^2 * (d.^(k-1) * transpose(d.^(k-1)) - k * (d.^(k-1) * transpose(d.^(k-1))) + (k-1)*diag(d.^(2*k)));
        eud_hessian = @(d,k) (k-1)*eud_func(d,k)/sum(d.^k)^2 * (-d.^(k-1) * transpose(d.^(k-1)) + sum(d.^k)*diag(d.^(k-2)));

switch method
    case 'int_gauss'
        expPowerSum = apm_powerSumProb(expDose,covDose,k,'int_gauss');
        [~,stdPowerSumTaylor] = apm_powerSumProb(expDose,covDose,k,'taylor_up');
        
        expEUD = real(1/fullVoxNum^(1/k) * apm_expPowerX(expPowerSum,stdPowerSumTaylor,1/k,'gauss'));
        raw2nd = real(1/fullVoxNum^(2/k) * apm_expPowerX(expPowerSum,stdPowerSumTaylor,2/k,'gauss'));

        stdEUD = sqrt(raw2nd - expEUD.^2);
        
    case 'taylor_steps'
        [expPowerSum,stdPowerSum] = apm_powerSumProb(expDose,covDose,k,'taylor_up');
        
        sqrt_func = @(d) 1 /fullVoxNum^(1/k) * d.^(1 / k);
        sqrt_grad = @(d) 1/fullVoxNum^(1/k) * 1/k * d.^(1/k - 1);
        sqrt_hess = @(d) 1/fullVoxNum^(2/k) * (k-1)/k^2 * d.^(1/k -2);
        
        [expEUD,stdEUD] = apm_taylorUncertaintyPropagation(expPowerSum,stdPowerSum,sqrt_func,sqrt_grad,sqrt_hess);
        %[expEUD,stdEUD] = apm_taylorUncertaintyPropagation(expPowerSum,stdPowerSum.^2,sqrt_func);
    case 'taylor_direct'
        %when using taylor expansion (single-step) 
        [expEUD,stdEUD] = apm_taylorUncertaintyPropagation(expDose,covDose,@(d) eud_func(d,k),@(d) eud_grad(d,k),@(d) eud_hessian(d,k));
    case 'taylor_direct_first'
        [expEUD,stdEUD] = apm_taylorUncertaintyPropagation(expDose,covDose,@(d) eud_func(d,k),@(d) eud_grad(d,k));%,@(d) eud_hessian(d,k)); %no 2nd order
        %[expEUD,stdEUD] = apm_taylorUncertaintyPropagation(expDose,covDose,@(d) eud_func(d,k)); %symbolic gradients
        
    case 'int_gamma'
        expPowerSum = apm_powerSumProb(expDose,covDose,k,'int_gamma');
        [~,stdPowerSumTaylor] = apm_powerSumProb(expDose,covDose,k,'taylor_up');

        expEUD = 1/fullVoxNum^(1/k) * apm_expPowerX(expPowerSum,stdPowerSumTaylor,1/k,'gamma');
        
        raw2nd = 1/fullVoxNum^(2/k) * apm_expPowerX(expPowerSum,stdPowerSumTaylor,2/k,'gamma');
        stdEUD = sqrt(raw2nd - expEUD.^2);
        
    case 'int_lognormal'
        [expPowerSum,stdPowerSum] = apm_powerSumProb(expDose,covDose,k,'int_lognormal');
        %[expPowerSum,stdPowerSum] = apm_powerSumProb(expDose,covDose,k,'int_lognormal');
        
        
        expEUD = 1/fullVoxNum^(1/k) * apm_expPowerX(expPowerSum,stdPowerSum,1/k,'lognormal');
        raw2nd = 1/fullVoxNum^(2/k) * apm_expPowerX(expPowerSum,stdPowerSum,2/k,'lognormal');

        stdEUD = sqrt(raw2nd - expEUD.^2);       
    otherwise
        error(['Method ''' method ''' not defined!']);
end
%expEUD = expEUD / nFac;
%stdEUD = stdEUD / nFac;
end

