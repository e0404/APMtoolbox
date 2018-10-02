function j = apm_cJacob(dij,w,vois)

d = dij*w;

DoseProjection          = sparse([]);

for v=1:numel(vois)
    doseVoi = d(vois(v).ix);
    for j = 1:numel(vois(v).cFunc)
        obj = vois(v).cFunc{j};
        if isa(obj,'DoseConstraints.matRad_DoseConstraint')
            cNew = [obj.computeDoseConstraintFunction(doseVoi) - obj.upperBounds(); obj.lowerBounds() - obj.computeDoseConstraintFunction(doseVoi)];
            valid = isfinite(cNew);
            jacobSub = obj.computeDoseConstraintJacobian(doseVoi);
            if valid(1)                
                DoseProjection  = [DoseProjection,sparse(find(vois(v).ix),1,jacobSub,numel(d),1)];
            end
            if valid(2)
                jacobSub = -jacobSub;
                DoseProjection  = [DoseProjection,sparse(find(vois(v).ix),1,jacobSub,numel(d),1)];
            end
        end
    end
end
if ~isempty(DoseProjection)
    j = transpose(DoseProjection' * dij);
else
    j = [];
end
%j = dij' * DoseProjection;
end

