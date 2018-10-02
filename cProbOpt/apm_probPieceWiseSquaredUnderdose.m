classdef apm_probPieceWiseSquaredUnderdose
    %apm_PROBPIECEWISESQUARED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties        
        dMin
        p
    end
    
    methods
        function obj = apm_probPieceWiseSquaredUnderdose(p,dMin)
           obj.dMin = dMin;
           obj.p = p;            
        end
        
        
        
        
        function f = objectiveFunction(obj,expDose,covDose)
            varDose = diag(covDose);
            stdDose = sqrt(varDose);
            
            muTransMin = obj.dMin - expDose;
            %muTransMin = obj.dMin - expDose;


            f = obj.p/numel(expDose) * sum(arrayfun(@(mu,sig) (mu^2 + sig^2) * obj.standardnormcdf(mu/sig) + mu*sig^2*obj.zeronormpdf(mu,sig),muTransMin,stdDose));
            %minSq = sum(p.*arrayfun(@(mu,sig) (mu^2 + sig^2) * standardnormcdf(mu/sig) + mu*sig^2*zeronormpdf(mu,sig),muTransMin,stdDose));
        end
        
        function fGradMu = gradientMeanDose(obj,expDose,covDose)
            varDose = diag(covDose);
            stdDose = sqrt(varDose);
            
            muTransMin = obj.dMin - expDose;
            
            %fGradMu = obj.p/numel(expDose) * arrayfun(@(mu,sig) 2*mu * obj.standardnormcdf(mu/sig) + mu*sig^2*obj.zeronormpdf(mu,sig),expDose,stdDose);           
            %fGradMu = obj.p/numel(expDose) * arrayfun(@(mu,sig) obj.zeronormpdf(mu,sig)*sig^2 + mu*erfc(-mu/(sqrt(2)*sig)),expDose,stdDose);
            fGradMu = -obj.p/numel(expDose) * arrayfun(@(mu,sig)  mu + 2*obj.zeronormpdf(mu,sig)*sig^2 + mu*erf(mu/sqrt(2*sig^2)),muTransMin,stdDose);
           
        end
        
        function fGradCov = gradientCovDose(obj,expDose,covDose)
            varDose = diag(covDose);
            fGradVar = obj.gradientVarDose(expDose,varDose);            
            fGradCov = diag(fGradVar);
        end
        
        function fGradVar = gradientVarDose(obj,expDose,varDose)
            muTransMin = obj.dMin - expDose;
            fGradVar = obj.p/numel(expDose) * arrayfun(@(mu,sig) obj.standardnormcdf(mu/sig) + mu^3 * obj.zeronormpdf(mu,sig) /(2*sig),muTransMin,sqrt(varDose));
        end
        
    end
    methods (Static)
        function val = standardnormcdf(x)
            val =  0.5*(1 + erf(x/sqrt(2)));
        end
        
        function val = zeronormpdf(mu,sig)
            val = 1/(sqrt(2*pi)*sig) * exp(-0.5 * (mu / sig)^2);
        end
    end
    
end

