function g = apm_probObjGrad(edij,cijlm,w,vois)

%muGrad = zeros(size(mu_d));
%covGrad = zeros(size(cov_d));

g = zeros(size(w));

mu_d = edij*w;
cov_d = apm_calcCovDose(cijlm,w);

for v=1:numel(vois)
    mu_d_voi = mu_d(vois(v).ix); %edij(vois(v).ix,:)*w;
    cov_d_voi = cov_d(vois(v).ix,vois(v).ix); %apm_calcCovDose(cijlm(vois(v).ix,:,vois(v).ix,:),w);
    for j = 1:numel(vois(v).probObjFunc)
        obj = vois(v).probObjFunc{j};
        
        gMu = edij(vois(v).ix,:)' * obj.gradientMeanDose(mu_d_voi,cov_d_voi);
        
        covGrad_voi = covGrad_voi(vois(v).ix,:,vois(v).ix);
        %covGrad_voi = apm_calcCovGrad(cijlm(vois(v).ix,:,vois(v).ix,:),w);
        
        objGradCov = obj.gradientCovDose(mu_d_voi,cov_d_voi);
        
        if isa(covGrad,'sptensor')
            %Product with tensor toolbox
            gCov = double(ttt(covGrad_voi,tensor(objGradCov),[1 3],[1 2])); 
        else
            %Product with reshape
            covGradMatrificated = reshape(permute(covGrad_voi,[2 1 3]),[size(covGrad_voi,2) size(covGrad_voi,1)^2]);
            objGradCovVectorized = reshape(objGradCov,[numel(objGradCov) 1]);
            gCov = covGradMatrificated * objGradCovVectorized;
        end
        
               
        
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

