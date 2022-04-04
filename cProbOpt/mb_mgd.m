function r = mb_mgd(cases,mu,sigma)

[m1 n1] = size(mu);
c = max([m1 n1]);
if m1 .* n1 ~= c
   error('Mu must be a vector.');
end

[m n] = size(sigma);
if m ~= n
   error('Sigma must be square');
end

if m ~= c
   error('The length of mu must equal the number of rows in sigma.');
end

[T p] = chol(sigma);
if p ~= 0
  error('Sigma must be a positive definite matrix.');
end

 
if m1 == c
  mu = mu';
end

mu = mu(ones(cases,1),:);

r = randn(cases,c) * T + mu;

end
