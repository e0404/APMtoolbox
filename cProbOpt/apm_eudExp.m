function expEUD = apm_eudExp(expDose,stdDose,k,method)
switch method
    case 'int_gauss'
        expPowerSum = apm_powerSumProb(expDose,covDose,k,'int_gauss');
        [~,stdPowerSumTaylor] = apm_powerSumProb(expDose,covDose,k,'taylor_up');
        
        expEUD = real(1/numel(expDose)^(1/k) * apm_expPowerX(expPowerSum,stdPowerSumTaylor,1/k,'gauss'));
        raw2nd = real(1/numel(expDose)^(2/k) * apm_expPowerX(expPowerSum,stdPowerSumTaylor,2/k,'gauss'));

        stdEUD = sqrt(raw2nd - expEUD.^2);
        
    case 'taylor_steps'
        [expPowerSum,stdPowerSum] = apm_powerSumProb(expDose,covDose,k,'taylor_up');
        
        sqrt_func = @(d) 1 / numel(expDose)^(1/k) * d.^(1 / k);
        sqrt_grad = @(d) 1/numel(expDose)^(1/k) * 1/k * d.^(1/k - 1);
        sqrt_hess = @(d) 1/numel(expDose)^(2/k) * (k-1)/k^2 * d.^(1/k -2);
        
        %[expEUD,stdEUD] = apm_taylorUncertaintyPropagation(expPowerSum,stdPowerSum,sqrt_func,sqrt_grad,sqrt_hess);
        [expEUD,stdEUD] = apm_taylorUncertaintyPropagation(expPowerSum,stdPowerSum.^2,sqrt_func);
    case 'taylor_direct'
        %when using taylor expansion (single-step) 

        eud_func = @(d,k) (sum(d.^k)/numel(d))^(1/k);
        eud_grad = @(d,k) transpose(eud_func(d,k)/sum(d.^k) * d.^(k-1));
        eud_hessian = @(d,k) eud_func(d,k)/sum(d.^k)^2 * (d.^(k-1) * transpose(d.^(k-1)) - k * (d.^(k-1) * transpose(d.^(k-1))) + (k-1)*diag(d.^(2*k)));
        
        [expEUD,stdEUD] = apm_taylorUncertaintyPropagation(expDose,covDose,@(d) eud_func(d,k),@(d) eud_grad(d,k));%,@(d) eud_hessian(d,k));
        %[expEUD,stdEUD] = apm_taylorUncertaintyPropagation(expDose,covDose,@(d) eud_func(d,k));%,@(d) eud_grad(d,k),@(d) eud_hessian(d,k));
        
    case 'int_gamma'
        expPowerSum = apm_powerSumProb(expDose,covDose,k,'int_gamma');
        [~,stdPowerSumTaylor] = apm_powerSumProb(expDose,covDose,k,'taylor_up');

        expEUD = 1/numel(expDose)^(1/k) * apm_expPowerX(expPowerSum,stdPowerSumTaylor,1/k,'gamma');
        
        raw2nd = 1/numel(expDose)^(2/k) * apm_expPowerX(expPowerSum,stdPowerSumTaylor,2/k,'gamma');
        stdEUD = sqrt(raw2nd - expEUD.^2);
        
    case 'int_lognormal'
        [expPowerSum,stdPowerSum] = apm_powerSumProb(expDose,covDose,k,'int_lognormal');
        %[~,stdPowerSumTaylor] = apm_powerSumProb(expDose,covDose,k,'taylor_up');
        
        expEUD = 1/numel(expDose)^(1/k) * apm_expPowerX(expPowerSum,stdPowerSum,1/k,'lognormal');
        raw2nd = 1/numel(expDose)^(2/k) * apm_expPowerX(expPowerSum,stdPowerSum,2/k,'lognormal');

        stdEUD = sqrt(raw2nd - expEUD.^2);
    otherwise
        error(['Method ''' method ''' not defined!']);
end


end

