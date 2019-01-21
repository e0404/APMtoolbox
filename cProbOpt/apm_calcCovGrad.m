function covGrad = apm_calcCovGrad(Vijlm,w)

if isa(Vijlm,'sptensor')
    covGrad = 2 * ttv(Vijlm,w,4);
else
    covGrad = 2 * permute(reshape(reshape(permute(Vijlm,[1 3 2 4]),size(Vijlm,1)^2*numel(w),numel(w)) * w,size(Vijlm,1),size(Vijlm,3),numel(w)),[1 3 2]);
end 

end

