eud_func = @(d,k) (sum(d.^k)/numel(d))^(1/k);
eud_grad = @(d,k) transpose(eud_func(d,k)/sum(d.^k) * d.^(k-1));
eud_hessian = @(d,k) (k-1)*eud_func(d,k)/sum(d.^k)^2 * (-d.^(k-1) * transpose(d.^(k-1)) + sum(d.^k)*diag(d.^(k-2)));
%eud_taylorApprox = @(x,k,x0) eud_func(x0,k) + eud_grad(x0,k)*(x-x0) + 0.5*(x-x0)' * eud_hessian(x0,k) * (x-x0);
eud_taylorApprox = @(x,k,x0) eud_func(x0,k) + eud_grad(x0,k)*(x-x0) + 0.5*(x-x0)' * eud_hessian(x0,k) * (x-x0);

nVar = 4;
nVals = 100;

k = 3.5;

xVals = rand(nVar,1)*10 + 1;

figure;
p=numSubplots(nVar);

for i = 1:nVar
    subplot(p(1),p(2),i);
    xSpace = linspace(xVals(i)-1,xVals(i)+1,nVals);    
    for step=1:nVals
        x_tmp = xVals;
        x_tmp(i) = xSpace(step);
        eudVals(step) = eud_func(x_tmp,k);
        eudVals_t(step) = eud_taylorApprox(x_tmp,k,xVals);
    end
    plot(xSpace,eudVals); hold on;
    plot(xSpace,eudVals_t);
    plot(xVals(i),eud_func(xVals,k),'o');
end
        


