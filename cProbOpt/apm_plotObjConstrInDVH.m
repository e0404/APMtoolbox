function [outputArg1,outputArg2] = apm_plotObjConstrInDVH(h,vois,prob)
if ~isvalid(h)
    warning('No valid axis provided! Skipping...');
end

hold(h,'on');

for v = 1:numel(vois)
    voi = vois(v);
    %Check if Probabilistic or Nominal
    if prob
        plotNominalObjectives(h,voi.probObjFunc);
        plotProbabilisticConstraints(h,voi.probCFunc);       
    else
        %Objectives
        plotNominalObjectives(h,voi.objFunc);   
        plotNominalConstraints(h,voi.cFunc);
    end
end
end

function plotNominalObjectives(h,objCell)
for o = 1:numel(objCell)
    obj = objCell{o};
    objName = class(obj);
    
    switch objName
        case 'DoseObjectives.matRad_SquaredDeviation'
            
        case 'DoseObjectives.matRad_SquaredUnderdosing'
            
        case 'DoseObjectives.matRad_SquaredOverdosing'
            
        case 'DoseObjectives.matRad_MinDVH'
            
        case 'DoseObjectives.matRad_MaxDVH'
            
        case 'DoseObjectives.matRad_MeanDose'
            
        case 'DoseObjectives.matRad_EUD'
            
        otherwise
    end
end
end

function plotNominalConstraints(h,constrCell)
for o = 1:numel(constrCell)
    constr = constrCell{o};
    constrName = class(constr);
    
    switch constrName
        case 'DoseConstraints.matRad_MinMaxDose'
            
        case 'DoseConstraints.matRad_MinMaxDVH'
            dP = constr.parameters{1};
            Vmin = constr.parameters{2};
            Vmax = constr.parameters{3};
            if Vmin > 0
                plot(h,dP,Vmin,'^k','MarkerFaceColor','k');            
            else
                plot(h,dP,Vmax,'vk','MarkerFaceColor','k');            
            end
        case 'DoseConstraints.matRad_MinMaxEUD'
            
        case 'DoseConstraints.matRad_MinMaxMeanDose'

        otherwise
    end
end
end

function plotProbabilisticObjectives(h,probObjCell)

end

function plotProbabilisticConstraints(h,probConstrCell)
for o = 1:numel(probConstrCell)
    constr = probConstrCell{o};
    constrName = class(constr);
    
    switch constrName          
        case 'apm_probDvhMinQuantileConstraint'
            dP = constr.dParam;
            Vmin = constr.minV;
            p = constr.p;            
            plot(h,dP,Vmin,'^k','MarkerFaceColor','k'); 
        case 'apm_probDvhMaxQuantileConstraint'
            dP = constr.dParam;
            Vmax = constr.maxV;
            p = constr.p;            
            plot(h,dP,Vmax,'vk','MarkerFaceColor','k'); 
       
        otherwise
    end
end
end