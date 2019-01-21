function vois = apm_setDefaultObjectivesAndConstraints(vois,obj,constr)

if nargin < 3
    constr = '';
end

if nargin < 2
    obj = 'sqDev';
end

%obj = 'sqDev', 'pwSqDev' or 'pwSqDevFake'

dRef = max([vois(:).dPres]);

%set objective function
switch obj
    case 'sqDev'        
        %objFunc = @(x) (dose_ij*x - dPres)'*diag(p)*(dose_ij*x - dPres);
        %gradFunc = @(x) ((2*p.*(dose_ij*x - dPres))' * dose_ij)';
        for v=1:numel(vois)
            optFunc = DoseObjectives.matRad_SquaredDeviation;
            optFunc.parameters{1} = vois(v).dPres;
            optFunc.penalty = vois(v).p;
            vois(v).objFunc{end+1} = optFunc;        
            
            probOptFunc = apm_probLeastSquares(vois(v).p,vois(v).dPres); 
            vois(v).probObjFunc{end+1} = probOptFunc;  
        end      
    case 'pwSqDev'           
        for v=1:numel(vois)
            
            if vois(v).dObjMin == vois(v).dObjMax && vois(v).dObjMin > 0
                optFunc = DoseObjectives.matRad_SquaredDeviation;
                optFunc.parameters{1} = vois(v).dPres;
                optFunc.penalty = vois(v).p;
                vois(v).objFunc{end+1} = optFunc; 
                
                probOptFunc = apm_probLeastSquares(vois(v).p,vois(v).dPres); 
                vois(v).probObjFunc{end+1} = probOptFunc;  
            else
                if vois(v).dObjMin > 0
                    optFunc = DoseObjectives.matRad_SquaredUnderdosing;
                    optFunc.parameters{1} = vois(v).dObjMin;
                    optFunc.penalty = vois(v).p;
                    vois(v).objFunc{end+1} = optFunc;
                    
                    probOptFunc = apm_probPieceWiseSquaredUnderdose(vois(v).p,vois(v).dObjMin); 
                    vois(v).probObjFunc{end+1} = probOptFunc; 
                end
                optFunc = DoseObjectives.matRad_SquaredOverdosing;
                optFunc.parameters{1} = vois(v).dObjMax;
                optFunc.penalty = vois(v).p;
                vois(v).objFunc{end+1} = optFunc;
                
                probOptFunc = apm_probPieceWiseSquaredOverdose(vois(v).p,vois(v).dObjMax);
                vois(v).probObjFunc{end+1} = probOptFunc;
            end            
        end
    case 'pwSqDevFake'           
        for v=1:numel(vois)
            
            if vois(v).dObjMin == vois(v).dObjMax && vois(v).dObjMin > 0
                optFunc = DoseObjectives.matRad_SquaredDeviation;
                optFunc.parameters{1} = vois(v).dPres;
                optFunc.penalty = vois(v).p;
                vois(v).objFunc{end+1} = optFunc; 
                
                probOptFunc = apm_probLeastSquares(vois(v).p,vois(v).dPres); 
                vois(v).probObjFunc{end+1} = probOptFunc;  
            else
                if vois(v).dObjMin > 0
                    optFunc = DoseObjectives.matRad_SquaredUnderdosing;
                    optFunc.parameters{1} = vois(v).dObjMin;
                    optFunc.penalty = vois(v).p;
                    vois(v).objFunc{end+1} = optFunc;
                    
                    probOptFunc = apm_probFakePieceWiseSquaredUnderdose(vois(v).p,vois(v).dObjMin); 
                    vois(v).probObjFunc{end+1} = probOptFunc; 
                end
                optFunc = DoseObjectives.matRad_SquaredOverdosing;
                optFunc.parameters{1} = vois(v).dObjMax;
                optFunc.penalty = vois(v).p;
                vois(v).objFunc{end+1} = optFunc;
                
                probOptFunc = apm_probFakePieceWiseSquaredOverdose(vois(v).p,vois(v).dObjMax);
                vois(v).probObjFunc{end+1} = probOptFunc;
            end            
        end
    otherwise
        error(['Objective ''' obj ''' not implemented!']);
end


%constr = 'meanMax', 'meanMin', 'EUDmax', 'EUDmin', 'DVHmin', 'DVHmax',
%'maxDose', 'minDose'

switch constr
    case 'DVHmin'
        dvhMinVol = 0.98;
        dvhDparam = 0.98*dRef;  
        optFunc = DoseConstraints.matRad_MinMaxDVH;
        optFunc.parameters{1} = dvhDparam;
        optFunc.parameters{2} = dvhMinVol;
        optFunc.parameters{3} = Inf;
        optFunc.referenceScalingVal = 1e-9;
        %optFunc.voxelScalingRatio = 10;
        vois(1).cFunc{end+1} = optFunc;        
        %constrFunc = @(x) dvhMinVol - vois(1).optFunc{1}.computeDoseConstraintFunction(dose_ij(vois(1).ix,:)*x);
        %constrJacob = @(x) -transpose(transpose(matRad_jacobFunc(dose_ij(vois(1).ix,:)*x,constraint,dvhDparam)) * dose_ij(vois(1).ix,:));
        %constrJacob = @(x) -transpose(transpose(vois(1).optFunc{1}.computeDoseConstraintJacobian(dose_ij(vois(1).ix,:)*x)) * dose_ij(vois(1).ix,:));
        
        dvhMinProbability = 0.15865; %1 standard deviation
        cObj = apm_probDvhMinQuantileConstraint(dvhDparam,dvhMinVol,dvhMinProbability,'normal');
        vois(1).probCFunc{end+1} = cObj;
        
    case 'EUDmin'
        eudMin = 0.99*dRef ;
        %eudK = -20;
        vois(1).eudK;
        
        optFunc = DoseConstraints.matRad_MinMaxEUD;
        optFunc.parameters{1} = eudK;
        optFunc.parameters{2} = eudMin;
        optFunc.parameters{3} = Inf;
        vois(1).cFunc{end+1} = optFunc;
        
        eudMinProbability = 0.15865; %1 standard deviation
        %eudMinProbability = 0.5;
        cObj = apm_probEudMinQuantileConstraint(eudMin,eudK,eudMinProbability,'normal');
        vois(1).probCFunc{end+1} = cObj;                
    case 'minDose'
        minDose = 0.95*vois(1).dPres;
        optFunc = DoseConstraints.matRad_MinMaxDose;
        optFunc.parameters{1} = minDose;
        optFunc.parameters{2} = Inf;        
        vois(1).cFunc{end+1} = optFunc;
        
        eudMinProbability = 0.05;
        cObj = apm_probEudMinQuantileConstraint(minDose,-100,eudMinProbability,'normal');
        vois(1).probCFunc{end+1} = cObj; 
    case 'DVHmax'
        dvhMaxVol = 0.3;
        dvhDparam = 0.5*vois(1).dPres;
        optFunc = DoseConstraints.matRad_MinMaxDVH;
        optFunc.parameters{1} = dvhDparam;
        optFunc.parameters{2} = -Inf;
        optFunc.parameters{3} = dvhMaxVol;
        optFunc.referenceScalingVal = 1e-9;
        vois(2).cFunc{end+1} = optFunc;   
        
        %dvhMaxProbability = 1-0.15865; %1 standard deviation
        dvhMaxProbability = 0.997; %1 standard deviation
        cObj = apm_probDvhMaxQuantileConstraint(dvhDparam,dvhMaxVol,dvhMaxProbability,'normal');
        vois(2).probCFunc{end+1} = cObj;
        
    case 'EUDmax'
        eudMax = 0.4*dRef ;
        eudK = vois(2).eudK;
        optFunc = DoseConstraints.matRad_MinMaxEUD;
        optFunc.parameters{1} = eudK;
        optFunc.parameters{2} = -Inf;
        optFunc.parameters{3} = eudMax;
        vois(2).cFunc{end+1} = optFunc;    
        
        %eudMaxProbability = 1-0.15865; %1 standard deviation
        %eudMaxProbability = 0.5; 
        eudMaxProbability = 0.95;
        cObj = apm_probEudMaxQuantileConstraint(eudMax,eudK,eudMaxProbability,'normal');
        vois(2).probCFunc{end+1} = cObj;
    case 'meanMax'
        eudMax = 0.4*dRef ;
        eudK = 1;
        optFunc = DoseConstraints.matRad_MinMaxEUD;
        optFunc.parameters{1} = eudK;
        optFunc.parameters{2} = -Inf;
        optFunc.parameters{3} = eudMax;
        vois(2).cFunc{end+1} = optFunc;    
        
        %eudMaxProbability = 1-0.15865; %1 standard deviation
        %eudMaxProbability = 0.5; 
        eudMaxProbability = 0.95;
        cObj = apm_probEudMaxQuantileConstraint(eudMax,eudK,eudMaxProbability,'normal');
        vois(2).probCFunc{end+1} = cObj;
    case 'maxDose'
        maxDose = 0.75*dRef ;
        optFunc = DoseConstraints.matRad_MinMaxDose;
        optFunc.parameters{1} = -Inf;
        optFunc.parameters{2} = maxDose;        
        vois(2).cFunc{end+1} = optFunc;       
        
        %maxDoseProbability = 1-0.15865;
        %maxDoseProbability = 1-0.5;
        maxDoseProbability = 0.95;
        cObj = apm_probEudMaxQuantileConstraint(maxDose,100,maxDoseProbability,'normal');
        vois(2).probCFunc{end+1} = cObj;
    otherwise       
end
end

