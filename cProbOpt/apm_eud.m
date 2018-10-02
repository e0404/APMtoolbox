function eud = apm_eud(d,exponent)
    eud = (1/numel(d) * sum(d.^exponent))^(1/exponent);
end

