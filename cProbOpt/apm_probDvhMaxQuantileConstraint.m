classdef apm_probDvhMaxQuantileConstraint
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dParam 
        maxV
        p
        dist
        considerCovGradSymmetry = false;
    end
    
    properties (Access = private)
        dvExpBuffered
        dvStdBuffered       
    end
    
    methods
        function obj = apm_probDvhMaxQuantileConstraint(dParam,maxVolume,probability,dist)           
            obj.dParam = dParam; 
            obj.maxV = maxVolume;
            
            if nargin < 4
                dist = 'normal';
            end
            obj.dist = dist;
            
            if nargin < 3
                probability = 0.5;
            end
            obj.p = probability;            
        end
        
        
        function c = constraintFunction(obj,expDose,covDose)
            [dvExp,dvStd] = apm_doseVolumeProb(expDose,covDose,obj.dParam);
            %[dvExp] = apm_doseVolumeProb(expDose,covDose,dParam);
            
            switch obj.dist
                case 'normal'
                    finv = dvExp + dvStd*sqrt(2)*erfinv(2*obj.p - 1);
                    %case 'gamma'
                    %    [p,b] = apm_transformMeanStdToGammaParameters(dvExp,dvStd);
            end
            
            c = finv-obj.maxV;
            
            %[expGrad,varGrad] = probGrads(expDose,covDose,dParam);
        end
        
        function muGrad = gradientMeanDose(obj,expDose,covDose)
            [dvExp,dvStd] = apm_doseVolumeProb(expDose,covDose,obj.dParam);
            
            varDose = diag(covDose);
            stdDose = sqrt(varDose);
            
            switch obj.dist
                case 'normal'                         
                    shift = (expDose-obj.dParam) ./ stdDose;
                    Nvals = obj.normpdf(obj.dParam,expDose,stdDose);
                    
                    %Gradient of the expectation value of DV
                    muGradExp = Nvals ./ numel(expDose);
                    
                    %probability term
                    %2nd raw moment derivative first (1-Fbvn)
                    iSum = zeros(size(expDose));
                    for i = 1:numel(expDose)                        
                        ix = logical(ones(size(expDose)));
                        ix(i) = false;
                        covRow = covDose(ix,i);
                        rhoRow = covRow ./(stdDose(i)*stdDose(ix));
                        iSum(i) = 1+2*sum(arrayfun(@(xl,rho) obj.normcdf(xl,rho*shift(i),sqrt(1-rho^2)),shift(ix),rhoRow));    
                    end
                    
                    %Now we have the gradient of the 2nd raw moment
                    rawMoment2gradMu = Nvals .* iSum;
                    rawMoment2gradMu = rawMoment2gradMu ./ numel(expDose)^2;
                    
                    %squared expectation
                    squaredExpGrad = 2 * dvExp * muGradExp; %Gradient has been verified
                                                
                    %Now we have the gradient of the 2nd central moment
                    muGradVar = rawMoment2gradMu - squaredExpGrad;  
                    
                    muGrad = muGradExp + 0.5*sqrt(2)*erfinv(2*obj.p - 1)/dvStd * muGradVar;
            end
        end
        
        function varGradExp = gradientVarDose(obj,expDose,varDose)
            switch obj.dist
                case 'normal'
                    normpdf = @(x,mu,sig) exp(-0.5*(x - mu).^2 ./ sig.^2)./(sqrt(2*pi) * sig);
                    
                    stdDose = sqrt(varDose);
                    
                    Nvals = normpdf(obj.dParam,expDose,stdDose);

                    varGradExp = (obj.dParam - expDose)./(2*stdDose.^2) .* Nvals ./ numel(expDose);
            end
        end
        
        
        function covGrad = gradientCovDose(obj,expDose,covDose)
            [dvExp,dvStd] = apm_doseVolumeProb(expDose,covDose,obj.dParam);
            
            varDose = diag(covDose);
            stdDose = sqrt(varDose);
                                
            covGrad = zeros(size(covDose));
            
            switch obj.dist
                case 'normal'
                    bivarPdf = @(x,mu,Sig) 1/(2*pi*sqrt(det(Sig))) * exp(-0.5*(x - mu)' / Sig * (x-mu));
                    
                    shift = -(obj.dParam - expDose) ./ stdDose;
                    
                    Nvals = obj.normpdf(obj.dParam,expDose,stdDose);
                    
                    %Prefactor for diagonals
                    varGrad = (-shift ./ (2*stdDose) .* Nvals) ./ numel(expDose);
                    
                    for i = 1:numel(expDose)
                        for l = i+1:numel(expDose)
                            %Off diagonal elements
                            covGrad(i,l) = 1/numel(covDose) * bivarPdf([obj.dParam; obj.dParam],[expDose(i); expDose(l)], [covDose(i,i), covDose(i,l); covDose(i,l), covDose(l,l)]);
                            covGrad(l,i) = covGrad(i,l);
                        end
                        
                        %Diagonal elements, sum of 1D cdfs
                        covRow = covDose(:,i);
                        rhoRow = covRow ./ (stdDose(i) * stdDose);
                        
                        
                        ix = logical(ones(size(expDose)));
                        ix(i) = false;
                        
                        covGrad(i,i) = varGrad(i)./numel(expDose) *(1+2*sum(arrayfun(@(xl,rho) obj.normcdf(xl,rho*shift(i),sqrt(1-rho^2)),shift(ix),rhoRow(ix))));                                                             
                    end                    
                    
                    %covGradVar = covGradVar ./ numel(covDose);
                    
                    squaredExpGrad = 2 * apm_doseVolumeProb(expDose,covDose,obj.dParam) * varGrad;
                    
                    %Compute gradient of the squared expectation value             
                    covGrad = covGrad - diag(squaredExpGrad);
                    
                    covGrad = diag(varGrad) + 0.5*sqrt(2)*erfinv(2*obj.p - 1)/dvStd * covGrad;
                    
                    %Take care of symmetry?
                    if obj.considerCovGradSymmetry
                        covGrad = 2*covGrad - diag(diag(covGrad));
                    end
                     
                    %covGradVar = -covGradVar;
            end
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
                x = (mu-x)./sig;
            end
            out = 0.5*erfc(x/sqrt(2));
        end        
    end
end

