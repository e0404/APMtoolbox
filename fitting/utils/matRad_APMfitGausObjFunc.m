function [ F, da, fa ] = matRad_APMfitGausObjFunc(a,x,y,returnRowVector)

%objective function

n = numel(a);
m = numel(y);

if mod(n,3) ~= 0
    error('inconsistent input\n');
end

if iscolumn(a)
    a = a';
end
w   = a(1:n/3);
mu  = a(n/3+1:2*n/3);
ell = a(2*n/3+1:end);

q = bsxfun(@minus,x,mu');
r = q.^2;
s = ell'*ones(1,m);
t = 1./sqrt(2*pi)./s .* exp(-r ./ (2*s.^2));


fa = w*t;
deviation = fa-y;
F = sum(deviation.^2);

dw   = 2 * t * deviation';
dmu  = 2 * ((w'*ones(1,m)) .* q ./ s.^2 .* t) * deviation';
dell = 2* ( w'*ones(1,m)) .* (-t ./ (ell'*ones(1,m)) + t .* r ./ (ell.^3'*ones(1,m))) * deviation';

if returnRowVector
    da = [dw' dmu' dell'];
else
    da = [dw' dmu' dell']';
    fa = fa';
end

end

