function c = apm_cFunc(dij,w,vois)

d = dij*w;

c = [];

for v=1:numel(vois)
    doseVoi = d(vois(v).ix);
    for j = 1:numel(vois(v).cFunc)
        obj = vois(v).cFunc{j};
        if isa(obj,'DoseConstraints.matRad_DoseConstraint')
            cNew = [obj.computeDoseConstraintFunction(doseVoi) - obj.upperBounds(); obj.lowerBounds() - obj.computeDoseConstraintFunction(doseVoi)];
                
            c = [c; cNew(isfinite(cNew))];            
        end
    end
end

end

