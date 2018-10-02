classdef apm_probFakePieceWiseSquaredUnderdose
    %apm_PROBPIECEWISESQUARED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties        
        dMax
        p
    end
    
    methods
        function obj = apm_probFakePieceWiseSquaredUnderdose(p,dMin)
           obj.dMin = dMin;
           obj.p = p;            
        end
        
        
        
        
        function f = objectiveFunction(obj,expDose,covDose)           
            dev = expDose - obj.dMin;
            %muTransMin = obj.dMin - expDose;
            
            %positivity
            dev(dev > 0) = 0;
           
            % claculate objective function
            f = obj.p/numel(expDose) * (trace(covDose) + (dev'*dev));
            
       
            
            %minSq = sum(p.*arrayfun(@(mu,sig) (mu^2 + sig^2) * standardnormcdf(mu/sig) + mu*sig^2*zeronormpdf(mu,sig),muTransMin,stdDose));
        end
        
        function fGradMu = gradientMeanDose(obj,expDose,covDose)           
            dev = expDose - obj.dMin;
            %muTransMin = obj.dMin - expDose;
            
            %positivity
            dev(dev > 0) = 0;
            
            % calculate delta
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

