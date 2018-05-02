function covGrad = apm_calcCovGrad(Vijlm,w)

%covGrad = permute(reshape(reshape(permute(Vijlm,[1 3 4 2]),size(Vijlm,1)^2*numel(w),numel(w)) * w,size(Vijlm,1),size(Vijlm,3),numel(w)),[1 3 2]);
%covGrad = etprod('ijl',Vijlm,'ijlm',w,'j');
covGrad = 2*double(ttv(tensor(Vijlm),w,4));


end

