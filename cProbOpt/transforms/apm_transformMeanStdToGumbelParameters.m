function [epsilon,alpha] = apm_transformMeanStdToGumbelParameters(mu, sigma)
% Return the parameters of a Gumbel distribution

% See Mahdi, Smail & Cenac, Myrtene. (2012). "Estimating Parameters of Gumbel Distribution using the Methods of Moments, probability weighted Moments and maximum likelihood". Revista de Matemática: Teoría y Aplicaciones. 12. 151. 10.15517/rmta.v12i1-2.259. 

J = 1.978;
gamma = 0.577215664901532;
alpha = sqrt(sigma.^2 ./ (J-gamma.^2));
epsilon = mu - gamma .* alpha;

end

