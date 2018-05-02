function eudHessian = apm_eudHessian(d,k)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


eudHessian = (k-1)*apm_eud(d,k)/sum(d.^k)^2 * (-d.^(k-1) * transpose(d.^(k-1)) + sum(d.^k)*diag(d.^(k-2)));

end

