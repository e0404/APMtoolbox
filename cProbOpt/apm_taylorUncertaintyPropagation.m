function [expT,stdT] = apm_taylorUncertaintyPropagation(mu,SIGMA,func,jacobFunc,hessianFunc)

onlyFirstOrder = false;
if nargin < 5
    onlyFirstOrder = true;
end

if nargin < 4
    %muSym = sym('muSym',size(mu));
    %f(muSym) = func(muSym);
    %fSym = feval(symengine, 'fp::unapply', f, P(1), P(2));
    %fG = jacobian(f);
    %fH = hessian(f);
    
    
    %jacobFunc = matlabFunction(fG);
    %hessianFunc = matlabFunction(fH);
    %onlyFirstOrder = false;
    
    hessEval = hessian(func,mu);
    gradEval = gradest(func,mu);
    onlyFirstOrder  = false;
else
    gradEval = jacobFunc(mu);
    if ~onlyFirstOrder
        hessEval = hessianFunc(mu);
    end
end

%hessianEval = hessianFunc(mu);

expT = func(mu);
if ~onlyFirstOrder
    %expT = expT + sum(sum(hessianFunc(mu)'*SIGMA*hessianFunc(mu)));
    order2 = 0.5 * sum(sum(SIGMA.*hessEval));
    %order2 = 0.5*trace(SIGMA*hessEval);
    expT = expT + order2;
end

stdT = sqrt(gradEval * SIGMA * gradEval');

end

