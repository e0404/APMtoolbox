function dv = apm_doseVolume(d,dParam)
    dv = numel(d(d >= dParam)) / numel(d);
end

