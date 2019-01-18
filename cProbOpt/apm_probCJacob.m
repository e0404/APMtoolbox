function jacob = apm_probCJacob(edij,cijlm,w,vois)

jacob = zeros(numel(w),0);

covGrad = apm_calcCovGrad(cijlm,w);
mu_d = edij*w;
cov_d = apm_calcCovDose(cijlm,w);

for v=1:numel(vois)
    mu_d_voi = mu_d(vois(v).ix); %edij(vois(v).ix,:)*w;
    cov_d_voi = cov_d(vois(v).ix,vois(v).ix); %apm_calcCovDose(cijlm(vois(v).ix,:,vois(v).ix,:),w);
    for j = 1:numel(vois(v).probCFunc)
        obj = vois(v).probCFunc{j};
        jvoi_mu = edij(vois(v).ix,:)' * obj.gradientMeanDose(mu_d_voi,cov_d_voi);
        
        
        objGradCov = obj.gradientCovDose(mu_d_voi,cov_d_voi);
        
        covGrad_voi = covGrad_voi(vois(v).ix,:,vois(v).ix);
        %covGrad_voi = apm_calcCovGrad(cijlm(vois(v).ix,:,vois(v).ix,:),w);
        
        if isa(covGrad,'sptensor')
            %Product with tensor toolbox
            jvoi_cov = double(ttt(covGrad_voi,tensor(objGradCov),[1 3],[1 2]));
        else
            %Product with reshape/permute
            
            covGradMatrificated = reshape(permute(covGrad_voi,[2 1 3]),[size(covGrad_voi,2) size(covGrad_voi,1)^2]); %Matrix BxV^2
            objGradCovVectorized = reshape(objGradCov,[numel(objGradCov) 1]); %vector V^2
            jvoi_cov = covGradMatrificated * objGradCovVectorized;
        end
        
        jacob(:,end+1) = jvoi_mu + jvoi_cov;
    end
end

end

% function covGrad = calcCovGrad(Vijlm,w)
%     covGrad = 2*double(ttv(tensor(Vijlm),w,4));
% end

