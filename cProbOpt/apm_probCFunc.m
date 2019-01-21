function c = apm_probCFunc(edij,cijlm,w,vois)

c = [];

mu_d = edij*w;
cov_d = apm_calcCovDose(cijlm,w);

for v=1:numel(vois)
    %mu_d_voi = edij(vois(v).ix,:)*w;
    %cov_d_voi = apm_calcCovDose(cijlm(vois(v).ix,:,vois(v).ix,:),w);
    mu_d_voi = mu_d(vois(v).ix);
    cov_d_voi = cov_d(vois(v).ix,vois(v).ix);
    
    for j = 1:numel(vois(v).probCFunc)
        obj = vois(v).probCFunc{j};
                    
        c = [c; obj.constraintFunction(mu_d_voi,cov_d_voi)];            
        
    end
end

end

