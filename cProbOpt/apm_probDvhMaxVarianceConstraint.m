classdef apm_probDvhMaxVarianceConstraint
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dParam 
        maxVariance
        dist
    end
    
    properties (Access = private)
        dvExpBuffered
        dvStdBuffered
        
        considerCovGradSymmetry = false;
    end
    
    methods
        function obj = apm_probDvhMaxVarianceConstraint(dParam,maxVariance,dist)           
            obj.dParam = dParam; 
            obj.maxVariance = maxVariance;
            
            if nargin < 3
                dist = 'normal';
            end
            obj.dist = dist;          
        end
        
        
        function c = constraintFunction(obj,expDose,covDose)
            [dvExp,dvStd] = apm_doseVolumeProb(expDose,covDose,obj.dParam);
            %[dvExp] = apm_doseVolumeProb(expDose,covDose,dParam);
                
            c = dvStd^2 - obj.maxVariance;                        
        end
        
        function muGradExp = gradientMeanDose(obj,expDose,covDose)
            %[~,dvStd] = apm_doseVolumeProb(expDose,covDose,obj.dParam);
            
            %cFuncWrapper = @(x) obj.constraintFunction(x,covDose);
            %muGradExp_est = transpose(gradest(cFuncWrapper,expDose));
            %return;
            
            switch obj.dist
                case 'normal'                         
                    varDose = diag(covDose);
                    stdDose = sqrt(varDose);
                    
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
                    %rawMoment2gradMu = muGradExp .* iSum;
                    %rawMoment2gradMu = muGradExp ./ numel(expDose) .* iSum;
                    rawMoment2gradMu = Nvals .* iSum;
                    rawMoment2gradMu = rawMoment2gradMu ./ numel(expDose)^2;
                    
                    %squared expectation
                    squaredExpGrad = 2 * apm_doseVolumeProb(expDose,covDose,obj.dParam) * muGradExp; %Gradient has been verified
                                                
                    %Now we have the gradient of the 2nd central moment
                    muGradExp = rawMoment2gradMu - squaredExpGrad;
                    %muGradExp = rawMoment2gradMu;
                    %muGradExp = squaredExpGrad;
                    %muGradExp = -muGradExp;                                                            
            end                   
        end
        
        function muGradExp = gradientMeanDose_varOnly(obj,expDose,varDose)
        	warning('Gradient needs Covariance Information, assuming no correlation!');
            muGradExp = obj.gradientMeanDose(obj,expDose,diag(varDose));            
        end
        
        function varGradVar = gradientVarDose(obj,expDose,varDose)
            warning('Gradient needs Covariance Information, assuming no correlation!');
            covGradVar = obj.gradientCovDose(obj,expDose,diag(varDose));
            varGradVar = diag(covGradVar);
        end
        
        function covGradVar = gradientCovDose(obj,expDose,covDose)
            %cFuncWrapper = @(x) obj.constraintFunction(expDose,x);
            %covGradVar = gradest(cFuncWrapper,covDose);
            %return;
            
            switch obj.dist
                case 'normal'
                    bivarPdf = @(x,mu,Sig) 1/(2*pi*sqrt(det(Sig))) * exp(-0.5*(x - mu)' / Sig * (x-mu));
                    
                    covGradVar = zeros(size(covDose));
                    
                    varDose = diag(covDose);
                    stdDose = sqrt(varDose);
                    
                    shift = -(obj.dParam - expDose) ./ stdDose;
                    
                    Nvals = obj.normpdf(obj.dParam,expDose,stdDose);
                    
                    %Prefactor for diagonals
                    varGrad = (-shift ./ (2*stdDose) .* Nvals) ./ numel(expDose);                    
                    
                    for i = 1:numel(expDose)
                        for l = i+1:numel(expDose)
                            %Off diagonal elements
                            covGradVar(i,l) = 1/numel(covDose) * bivarPdf([obj.dParam; obj.dParam],[expDose(i); expDose(l)], [covDose(i,i), covDose(i,l); covDose(i,l), covDose(l,l)]);
                            covGradVar(l,i) = covGradVar(i,l);
                        end
                        
                        %Diagonal elements, sum of 1D cdfs
                        covRow = covDose(:,i);
                        rhoRow = covRow ./ (stdDose(i) * stdDose);
                        
                        
                        ix = logical(ones(size(expDose)));
                        ix(i) = false;
                        
                        covGradVar(i,i) = varGrad(i)./numel(expDose) *(1+2*sum(arrayfun(@(xl,rho) obj.normcdf(xl,rho*shift(i),sqrt(1-rho^2)),shift(ix),rhoRow(ix))));                                                             
                    end                    
                    
                    %covGradVar = covGradVar ./ numel(covDose);
                    
                    squaredExpGrad = 2 * apm_doseVolumeProb(expDose,covDose,obj.dParam) * varGrad;
                    
                    %Compute gradient of the squared expectation value             
                    covGradVar = covGradVar - diag(squaredExpGrad);
                    
                    %Take care of symmetry?
                    if obj.considerCovGradSymmetry
                        covGradVar = 2*covGradVar - diag(diag(covGradVar));
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



