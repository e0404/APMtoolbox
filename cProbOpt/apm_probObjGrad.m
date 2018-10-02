function g = apm_probObjGrad(edij,cijlm,w,vois)

%muGrad = zeros(size(mu_d));
%covGrad = zeros(size(cov_d));

g = zeros(size(w));

for v=1:numel(vois)
    mu_d_voi = edij(vois(v).ix,:)*w;
    cov_d_voi = apm_calcCovDose(cijlm(vois(v).ix,:,vois(v).ix,:),w);
    for j = 1:numel(vois(v).probObjFunc)
        obj = vois(v).probObjFunc{j};
        
        gMu = edij(vois(v).ix,:)' * obj.gradientMeanDose(mu_d_voi,cov_d_voi);
        covGrad = apm_calcCovGrad(cijlm(vois(v).ix,:,vois(v).ix,:),w);
        objGradCov = obj.gradientCovDose(mu_d_voi,cov_d_voi);
        
        %Product with reshape
        covGradMatrificated = reshape(permute(covGrad,[2 1 3]),[size(covGrad,2) size(covGrad,1)^2]);
        objGradCovVectorized = reshape(objGradCov,[numel(objGradCov) 1]);
        gCov = covGradMatrificated * objGradCovVectorized;
        
        %Product with tensor toolbox
        %gCov = double(ttt(tensor(covGrad),tensor(objGradCov),[1 3],[1 2])); 
        
        g = g + (gMu + gCov);
        %muGrad(vois(v).ix) = muGrad(vois(v).ix) + obj.gradientMeanDose(mu_d_voi,cov_d_voi);
        %covGrad(vois(v).ix,vois(v).ix) = covGrad(vois(v).ix,vois(v).ix) + obj.gradientCovDose(mu_d_voi,cov_d_voi);            
    end
end


% chain rule
%g = zeros(numel(w),1);

%gMu = (muGrad' * edij)';
%gCov = double(ttt(tensor(calcCovGrad(cijlm,w)),tensor(covGrad),[1 3],[1 2]));

%g = gMu + gCov;
end

