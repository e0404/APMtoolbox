function covInfluenceEle = apm_calcCovInfluenceElement(i,j,l,m,x,mu,sigma,C_R,C_S,expDose_ij,f)


dev =  x([i l]) - mu([j m]);
%dev = mu([j k]) - x(i);

tmp_CR = C_R([j m],[j m]);

%Systematic Term
A = diag([sigma(j)^2 sigma(m)^2]) + diag(diag(tmp_CR)) + C_S([j m],[j m]);
fac = 1/(2*pi*sqrt(det(A)));

vUncorr = fac * exp(-.5 * dev' / A * dev);

%Add Random Term
A = A + flipud(diag(diag(flipud(tmp_CR)))); %Add off diagonals
fac = 1/(2*pi*sqrt(det(A)));

vCorr = fac * exp(-.5 * dev' / A * dev);

covInfluenceEleSys = vUncorr - expDose_ij(i,j)*expDose_ij(l,m);
covInfluenceEleRand = vCorr - vUncorr;

covInfluenceEle = 1/f * covInfluenceEleRand + covInfluenceEleSys;


end

