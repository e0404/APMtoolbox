function c = apm_probCFunc(edij,cijlm,w,vois)

c = [];

for v=1:numel(vois)
    mu_d_voi = edij(vois(v).ix,:)*w;
    cov_d_voi = apm_calcCovDose(cijlm(vois(v).ix,:,vois(v).ix,:),w);
    for j = 1:numel(vois(v).probCFunc)
        obj = vois(v).probCFunc{j};
                    
        c = [c; obj.constraintFunction(mu_d_voi,cov_d_voi)];            
        
    end
end

end

