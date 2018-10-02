function jacob = apm_probCJacob(edij,cijlm,w,vois)

jacob = zeros(numel(w),0);

for v=1:numel(vois)
    mu_d_voi = edij(vois(v).ix,:)*w;
    cov_d_voi = apm_calcCovDose(cijlm(vois(v).ix,:,vois(v).ix,:),w);
    for j = 1:numel(vois(v).probCFunc)
        obj = vois(v).probCFunc{j};
        jvoi_mu = edij(vois(v).ix,:)' * obj.gradientMeanDose(mu_d_voi,cov_d_voi);
        
        covGrad = apm_calcCovGrad(cijlm(vois(v).ix,:,vois(v).ix,:),w);
        objGradCov = obj.gradientCovDose(mu_d_voi,cov_d_voi);
        
        %Product with reshape/permute
        covGradMatrificated = reshape(permute(covGrad,[2 1 3]),[size(covGrad,2) size(covGrad,1)^2]); %Matrix BxV^2
        objGradCovVectorized = reshape(objGradCov,[numel(objGradCov) 1]); %vector V^2
        jvoi_cov = covGradMatrificated * objGradCovVectorized;
        
        %Product with tensor toolbox
        %jvoi_cov = double(ttt(tensor(covGrad),tensor(objGradCov),[1 3],[1 2]));
        
        jacob(:,end+1) = jvoi_mu + jvoi_cov;
    end
end

end

% function covGrad = calcCovGrad(Vijlm,w)
%     covGrad = 2*double(ttv(tensor(Vijlm),w,4));
% end

