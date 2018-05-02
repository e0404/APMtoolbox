function f = apm_objFunc(dij,w,vois)

d = dij*w;

f = 0;

for v=1:numel(vois)
    doseVoi = d(vois(v).ix);
    for j = 1:numel(vois(v).objFunc)
        obj = vois(v).objFunc{j};
        if isa(obj,'DoseObjectives.matRad_DoseObjective')
            f = f + obj.computeDoseObjectiveFunction(doseVoi);
        end
    end
end

