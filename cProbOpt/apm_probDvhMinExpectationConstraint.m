classdef apm_probDvhMinExpectationConstraint
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dParam 
        minV
        dist
    end
    
    properties (Access = private)
        dvExpBuffered
        dvStdBuffered
    end
    
    methods
        function obj = apm_probDvhMinExpectationConstraint(dParam,minVolume,dist)           
            obj.dParam = dParam; 
            obj.minV = minVolume;
            
            if nargin < 3
                dist = 'normal';
            end
            obj.dist = dist;          
        end
        
        
        function c = constraintFunction(obj,expDose,covDose)
            [dvExp,~] = apm_doseVolumeProb(expDose,covDose,obj.dParam);
            %[dvExp] = apm_doseVolumeProb(expDose,covDose,dParam);
                
            c = obj.minV - dvExp;                        
        end
        
        function muGradExp = gradientMeanDose(obj,expDose,covDose)
            %[~,dvStd] = apm_doseVolumeProb(expDose,covDose,obj.dParam);
            varDose = diag(covDose);
            muGradExp = obj.gradientMeanDose_varOnly(obj,expDose,varDose);                      
        end
        
        function muGradExp = gradientMeanDose_varOnly(obj,expDose,varDose)
           switch obj.dist
                case 'normal'                         
                    stdDose = sqrt(varDose);                    
                    Nvals = obj.normpdf(obj.dParam,expDose,stdDose);                    
                    %muGradExp = sum(Nvals) / numel(expDose);
                    muGradExp = -Nvals ./ numel(expDose);                    
            end
        end
        
        function varGradExp = gradientVarDose(obj,expDose,varDose)
            switch obj.dist
                case 'normal'
                    normpdf = @(x,mu,sig) exp(-0.5*(x - mu).^2 ./ sig.^2)./(sqrt(2*pi) * sig);
                    
                    stdDose = sqrt(varDose);
                    
                    Nvals = normpdf(obj.dParam,expDose,stdDose);

                    varGradExp = -(obj.dParam - expDose)./(2*stdDose.^2) .* Nvals ./ numel(expDose);
            end
        end
        
        function covGradExp = gradientCovDose(obj,expDose,covDose)
            covGradExp = diag(obj.gradientVarDose(obj,expDose,diag(covDose)));
        end
    end
    
    methods (Static)
        function out = normpdf(x,mu,sig)
            f = 1;
            if nargin > 1
                f = 1./sig;
                x = (x-mu).*f;                
            end
            out = exp(-0.5*x.^2)./sqrt(2*pi) .* f;
        end
        
        function out = normcdf(x,mu,sig)
            if nargin > 1
                x = (x-mu)./sig;
            end
            out = 0.5*(1+erf(x/sqrt(2)));
        end        
    end
end



