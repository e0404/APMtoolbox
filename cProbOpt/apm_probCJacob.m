function jacob = apm_probCJacob(edij,cijlm,w,vois)

jacob = zeros(numel(w),0);

for v=1:numel(vois)
    mu_d_voi = edij(vois(v).ix,:)*w;
    cov_d_voi = apm_calcCovDose(cijlm(vois(v).ix,:,vois(v).ix,:),w);
    for j = 1:numel(vois(v).probCFunc)
        obj = vois(v).probCFunc{j};
        jvoi_mu = edij(vois(v).ix,:)' * obj.gradientMeanDose(mu_d_voi,cov_d_voi);
        jvoi_cov = double(ttt(tensor(calcCovGrad(cijlm(vois(v).ix,:,vois(v).ix,:),w)),tensor(obj.gradientCovDose(mu_d_voi,cov_d_voi)),[1 3],[1 2]));
        
        jacob(:,end+1) = jvoi_mu + jvoi_cov;
    end
end

end

function covGrad = calcCovGrad(Vijlm,w)
    covGrad = 2*double(ttv(tensor(Vijlm),w,4));
end

