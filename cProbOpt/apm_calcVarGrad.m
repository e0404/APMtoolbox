function varGrad = apm_calcVarGrad(Vijm,w)
varGrad =  2*double(ttv(tensor(Vijm),w,3));
%varGrad =  2*reshape(reshape(permute(Vijm,[1 3 2]),size(Vijm,1)*numel(w),numel(w)) * w,size(Vijm,1),numel(w));
end

