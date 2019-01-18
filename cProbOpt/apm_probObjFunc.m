function f = apm_probObjFunc(edij,cijlm,w,vois)

mu_d = edij*w;
cov_d = apm_calcCovDose(cijlm,w);

f = 0;

for v=1:numel(vois)
    mu_d_voi = mu_d(vois(v).ix); %edij(vois(v).ix,:)*w;
    cov_d_voi = cov_d(vois(v).ix,vois(v).ix); %apm_calcCovDose(cijlm(vois(v).ix,:,vois(v).ix,:),w);
    for j = 1:numel(vois(v).probObjFunc)
        obj = vois(v).probObjFunc{j};
        f = f + obj.objectiveFunction(mu_d_voi,cov_d_voi);
    end
end