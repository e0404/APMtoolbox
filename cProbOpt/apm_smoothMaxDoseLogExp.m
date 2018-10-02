function smoothmax = apm_smoothMaxDoseLogExp(dose)
    smoothmax = log(sum(exp(dose)) - numel(dose) + 1);
end

