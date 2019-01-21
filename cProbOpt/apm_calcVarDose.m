function varDose = apm_calcVarDose(Vijm,w)

if isa(Vijm,'sptensor')
    varDose = ttt(Vijm,tensor(w*w'),[3 2],[1 2]);
else
    varDose = reshape(reshape(permute(Vijm,[1 3 2]),size(Vijm,1),numel(w)^2) * reshape(w*w',numel(w)^2,[]),size(Vijm,1),1); %Manual product via matricication
end

end

