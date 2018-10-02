function g = apm_objGrad(dij,w,vois)

d = dij*w;

dGrad = zeros(numel(d),1);

for v=1:numel(vois)
    doseVoi = d(vois(v).ix);
    for j = 1:numel(vois(v).objFunc)
        obj = vois(v).objFunc{j};
        if isa(obj,'DoseObjectives.matRad_DoseObjective')
            dGrad(vois(v).ix) = dGrad(vois(v).ix) + obj.computeDoseObjectiveGradient(doseVoi);
        end
    end
end


% chain rule
%g = zeros(numel(w),1);

g = (dGrad' * dij)';
end

