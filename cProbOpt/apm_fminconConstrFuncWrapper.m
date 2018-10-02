function [c,ceq,jc,jceq] = apm_fminconConstrFuncWrapper(x,cFunc,ceqFunc,cJacob,ceqJacob)

if isa(cFunc,'function_handle')
    c = cFunc(x);
else
    c = [];
end

if nargout > 1
    if isa(ceqFunc,'function_handle')
        ceq = ceqFunc(x);
    else
        ceq = [];
    end
end
if nargout > 2
    if isa(cJacob,'function_handle')
        jc = cJacob(x);
    else
        jc = [];
    end
end

if nargout > 3
    if isa(ceqJacob,'function_handle')
        jceq = ceqJacob(x);
    else
        jceq = [];
    end
end
end

