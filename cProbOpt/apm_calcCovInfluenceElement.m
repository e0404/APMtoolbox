function [covInfluenceEleSys , covInfluenceEleRand] = apm_calcCovInfluenceElement(i,j,l,m,x,spots,ucm,expDose_ij)



mu2D = [spots([j m]).mu]';
sigma2D = [spots([j m]).sigma]';

dev =  x([i l]) - mu2D;

%if false %any(abs(dev) > 5*sqrt(transpose(sigma2D).^2 + diag(ucm.covSys([j m],[j m])) + diag(ucm.covRand([j m],[j m]))))
%    covInfluenceEleSys = 0;
%    covInfluenceEleRand = 0;
%else


%dev = mu([j k]) - x(i);

tmp_CR = ucm.covRand([j m],[j m]);

%Systematic Term
A = diag(sigma2D.^2) + diag(diag(tmp_CR)) + ucm.covSys([j m],[j m]);
fac = 1/(2*pi*sqrt(det(A)));

vUncorr = fac * exp(-.5 * dev' / A * dev);

%Add Random Term
A = A + flipud(diag(diag(flipud(tmp_CR)))); %Add off diagonals
fac = 1/(2*pi*sqrt(det(A)));

vCorr = fac * exp(-.5 * dev' / A * dev);

covInfluenceEleSys = vUncorr - expDose_ij(i,j)*expDose_ij(l,m);
covInfluenceEleRand = vCorr - vUncorr;
%end
end

