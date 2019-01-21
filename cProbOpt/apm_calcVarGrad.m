function varGrad = apm_calcVarGrad(Vijm,w)

if isa(Vijm,'sptensor')
    varGrad =  2 * ttv(Vijm,w,3);
else
    varGrad =  2*reshape(reshape(permute(Vijm,[1 3 2]),size(Vijm,1)*numel(w),numel(w)) * w,size(Vijm,1),numel(w));
end

end

