function [f,gradF] = apm_fminconObjFuncWrapper(x,objFunc,gradFunc)
    f = objFunc(x);
    gradF = gradFunc(x);
end

