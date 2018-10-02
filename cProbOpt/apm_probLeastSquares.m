classdef apm_probLeastSquares
    %apm_PROBPIECEWISESQUARED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties        
        dRef
        p
    end
    
    methods
        function obj = apm_probLeastSquares(p,dRef)
           obj.dRef = dRef;
           obj.p = p;            
        end
        
        
        
        
        function f = objectiveFunction(obj,expDose,covDose)           
            dev = expDose - obj.dRef;
            %muTransMin = obj.dMin - expDose;

            %pDiag = obj.p/numel(expDose) * eye(size(covDose));
            
            %f = dev' * pDiag * dev + trace(pDiag*covDose);
            f = obj.p/numel(expDose) * (trace(covDose) + dev'*dev);
            %minSq = sum(p.*arrayfun(@(mu,sig) (mu^2 + sig^2) * standardnormcdf(mu/sig) + mu*sig^2*zeronormpdf(mu,sig),muTransMin,stdDose));
        end
        
        function fGradMu = gradientMeanDose(obj,expDose,covDose)           
            dev = expDose - obj.dRef;
            fGradMu = 2 * obj.p/numel(expDose) * dev;
           
        end
        
        function fGradCov = gradientCovDose(obj,expDose,covDose)
            fGradCov = obj.p/numel(expDose) * eye(size(covDose));
        end
        
        function fGradVar = gradientVarDose(obj,expDose,varDose)
            fGradVar = obj.p/numel(expDose) * ones(size(expDose));
        end
        
    end
    
end

