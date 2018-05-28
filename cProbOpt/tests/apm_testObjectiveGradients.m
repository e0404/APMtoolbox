%{
dParam = 0.95;
maxVar = 0.05^2;
n=4;

%cObj = apm_probDvhMaxVarianceConstraint(dParam,maxVar,'normal');
cObj = apm_probDvhMaxQuantileConstraint(dParam,0.3,0.25,'normal');
cObj.considerCovGradSymmetry = true;
%}
nVars=4;

k = -10;
%cObj = apm_probEudMaxQuantileConstraint(1,k,0.5);
fObj = apm_probPieceWiseSquaredOverdose(1,0.5);
%fObj = apm_probPieceWiseSquaredUnderdose(1,0.7);
%fObj = apm_probDvh

expDose = rand(nVars,1);
sigDose = rand(nVars,1);

nSteps = 100;

expDoseSpace = linspace(0,2,nSteps);

%covDose = eye(n);
%covDose = diag(sigDose.^2);
covDose = apm_createCovarianceMatrix(nVars,sigDose);
varDose = diag(covDose);
corrDose = covDose ./ (sigDose*sigDose');


%muGrad = cObj.gradientMeanDose(expDose,covDose);

varDoseSpaceRel = linspace(1e-6,3,nSteps);
corrSpaceRel = linspace(1e-6,1-1e-6,nSteps);


ix = 25;


hf = figure;
p = numSubplots(nVars);
for i=1:nVars
    subplot(p(1),p(2),i);
    expDose_tmp = expDose;
    
    for step=1:numel(expDoseSpace)
        expDose_tmp(i) = expDoseSpace(step);
        fSteps(step) = fObj.objectiveFunction(expDose_tmp,covDose);
        muGrad_tmp = fObj.gradientMeanDose(expDose_tmp,covDose);
        
        %cFuncWrapper = @(x) cObj.constraintFunction(x,covDose);
        fWrapper = @(x) fObj.objectiveFunction(x,covDose);
        muGrad_tmp_est = transpose(gradest(fWrapper,expDose_tmp));
        
        fDeriv(step) = muGrad_tmp(i);
        fDerivEst(step) = muGrad_tmp_est(i);
    end
    
    fNumGrad = gradient(fSteps,expDoseSpace);
    
    
    %plot(expDoseSpace,fSteps);hold on;
    %plot(expDoseSpace(ix),fSteps(ix),'ko');
    %tangent = fDeriv(ix) * (expDoseSpace - expDoseSpace(ix)) + fSteps(ix);
    %tangent_est = fDerivEst(ix) * (expDoseSpace - expDoseSpace(ix)) + fSteps(ix);
    %plot(expDoseSpace,tangent);
    %plot(expDoseSpace,tangent_est,'--');
    %xlim([min(expDoseSpace) max(expDoseSpace)]);
    plot(expDoseSpace,fNumGrad); hold on;
    plot(expDoseSpace,fDeriv);
    
    drawnow();
    %waitforbuttonpress;
end


hf = figure;

for i = 1:nVars    
    for l = 1:nVars
        subplot(nVars,nVars,nVars*(i-1)+l);
        covDose_tmp = covDose;
        varDose_tmp = varDose;
        
        for step = 1:nSteps            
            if i ~= l
                space_tmp = sigDose(i)*sigDose(l) * corrSpaceRel;
                covDose_tmp(i,l) = space_tmp(step);
                covDose_tmp(l,i) = covDose_tmp(i,l);
            else
                space_tmp = varDose(i)* varDoseSpaceRel;
                varDose_tmp(i) = space_tmp(step);
                covDose_tmp = (sqrt(varDose_tmp)*sqrt(varDose_tmp')) .* corrDose;
            end
            fSteps(step) = fObj.objectiveFunction(expDose,covDose_tmp);
            covGrad_tmp  = fObj.gradientCovDose(expDose,covDose_tmp);
            
            fDeriv(step) = covGrad_tmp(i,l);
        end
        
        %plot(space_tmp,fSteps);hold on;
        %plot(space_tmp(ix),fSteps(ix),'ko');
        %tangent = fDeriv(ix) * (space_tmp - space_tmp(ix)) + fSteps(ix);
        %tangent_est = cDerivEst(ix) * (space_tmp - space_tmp(ix)) + cSteps(ix);
        %plot(space_tmp,tangent);
        fNumGrad = gradient(fSteps,space_tmp);
        plot(space_tmp,fNumGrad); hold on;
        plot(space_tmp,fDeriv);
        
        
        %plot(space_tmp,tangent_est,'--');
        drawnow();
        %waitforbuttonpress;
    end
end

