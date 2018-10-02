function smoothmax = apm_smoothMaxDoseExpSum(dose,a)

%
smoothmax = sum(dose.*exp(a*dose)) / sum(exp(a*dose));


end

