function eudGrad = apm_eudGrad(d,k)

eudGrad = transpose(apm_eud(d,k)/sum(d.^k) * d.^(k-1));

end

