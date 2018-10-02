classdef apm_probEudMaxQuantileConstraint
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        k 
        maxEUD
        p
        dist
        considerCovGradSymmetry = false;
    end
    
    properties (Access = private)
        dvExpBuffered
        dvStdBuffered       
    end
    
    methods
        function obj = apm_probEudMaxQuantileConstraint(maxEUD,k,probability,dist)           
            obj.maxEUD = maxEUD; 
            obj.k = k;
            
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
            [eudExp,eudStd] = apm_eudProb(expDose,covDose,obj.k);            
            
            switch obj.dist
                case 'normal'
                    finv = eudExp + eudStd*sqrt(2)*erfinv(2*obj.p - 1);
            end            
            c = finv - obj.maxEUD;
        end
        
        function muGrad = gradientMeanDose(obj,expDose,covDose)
            %[eudExp,eudStd] = apm_eudProb(expDose,covDose,obj.k); 

            
            switch obj.dist
                case 'normal'                         
                    %exp. value term, First Order
                    muGrad = obj.eud_grad(expDose);
                    
                    %exp. value term, Second Order
                    
%                     hessGrad = zeros(size(muGrad));
%                     for i = 1:numel(expDose)
%                         hessDer = obj.eud_hessianMuGradDerivativeEle(expDose,i);
%                         hessGrad(i) = 0.5*sum(sum(covDose.*hessDer));
%                     end
%                     muGrad = muGrad + hessGrad;
                    
                    %muGrad = muGrad + diag(obj.eud_hessian(expDose));
                    % + diag(obj.eud_hessian(expDose));                                       
                    
            end
        end
        
        function varGrad = gradientVarDose(obj,expDose,varDose)
            covDose = diag(varDose);
            covGrad = obj.gradientCovDose(expDose,covDose);
            varGrad = diag(covGrad);
        end
        
        
        function covGrad = gradientCovDose(obj,expDose,covDose)
            [~,eudStd] = apm_eudProb(expDose,covDose,obj.k); 
            
            %varDose = diag(covDose);
            %stdDose = sqrt(varDose);
                                
            covGrad = zeros(size(covDose));
            
            switch obj.dist
                case 'normal'
                    covGradMu = 0.5*obj.eud_hessian(expDose);
                    
                    covGradVar = obj.eud_grad(expDose)*obj.eud_grad(expDose)';
                    
                    covGradAsym = covGradMu + 0.5*sqrt(2)*erfinv(2*obj.p - 1)/eudStd * covGradVar;
                    
                    if obj.considerCovGradSymmetry
                        for i=1:numel(expDose)
                            for l=1:numel(expDose)
                                structureMat = zeros(size(covGradAsym));
                                if i == l
                                    structureMat(i,:) = covDose(i,:)./(2*covDose(i,i));
                                    structureMat(:,i) = covDose(:,i)./(2*covDose(i,i));
                                end
                                structureMat(i,l) = 1;
                                structureMat(l,i) = 1;
                                
                                covGrad(i,l) = trace(covGradAsym*structureMat);                                
                            end
                        end
                        
                        %symmetric
                        %covGrad = 2*covGradAsym - diag(diag(covGradAsym));
                    else
                        covGrad = covGradAsym;
                    end
            end
        end
        
        function eud = eud(obj,d)
            eud = (1/numel(d) * sum(d.^obj.k))^(1/obj.k);
        end
        
        function gradEud = eud_grad(obj,d)
            gradEud = obj.eud(d)/sum(d.^obj.k) * d.^(obj.k-1);
        end
        
        function hessEud = eud_hessian(obj,d)
            hessEud = (obj.k-1)*obj.eud(d)/sum(d.^obj.k)^2 * (-d.^(obj.k-1) * transpose(d.^(obj.k-1)) + sum(d.^obj.k)*diag(d.^(obj.k-2)));
        end
        
        function der3EUD = eud_hessianMuGradDerivativeEle(obj,d,ix)
            toTheK = d.^obj.k;
            powersum = sum(toTheK);
            powersumMinusElement = powersum - toTheK;
            
            prefac = obj.eud(d)/powersum^3 * (obj.k-1) * d(ix)^(obj.k-2);
            
            der3EUD = prefac*ones(numel(d));
            
            for i = 1:numel(d)
                for l = i:numel(d)
                    if i == ix
                        if l ~= i
                            der3EUD(ix,l) = der3EUD(ix,l) * d(l)^(obj.k-1) * (d(ix)^obj.k * obj.k - powersumMinusElement(ix)*(obj.k-1));
                            der3EUD(l,ix) = der3EUD(ix,l);
                        else
                            der3EUD(ix,ix) = der3EUD(ix,ix) * powersumMinusElement(ix)/d(ix) * (powersumMinusElement(ix)*(obj.k-2) - d(ix)^obj.k*(obj.k+1));
                        end
                    else %We are on an element i,l ~= ix
                        if i == l %check if we are on a diagonal
                            der3EUD(i,l) = der3EUD(i,l) * d(ix) * d(i)^(obj.k-2) * (d(i)^obj.k * obj.k - powersumMinusElement(i)*(obj.k-1));
                        else
                            der3EUD(i,l) = der3EUD(i,l) * d(ix) * d(i)^(obj.k-1)*d(l)^(obj.k-1)*(2*obj.k-1);
                            der3EUD(l,i) = der3EUD(i,l);
                        end                        
                    end
                end
            end
        end
    end    
end

