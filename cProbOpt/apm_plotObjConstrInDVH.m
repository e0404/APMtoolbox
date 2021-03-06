function [outputArg1,outputArg2] = apm_plotObjConstrInDVH(h,vois,prob,plotBodyDVH)
if ~isvalid(h)
    warning('No valid axis provided! Skipping...');
end

if nargin < 4
    plotBodyDVH = true;
end

hold(h,'on');

for v = 1:numel(vois)
    voi = vois(v);
    if ~strcmp(vois(v).type,'BODY') || plotBodyDVH
        %Check if Probabilistic or Nominal
        if prob
            plotProbabilisticObjectives(h,voi.probObjFunc,voi.dvhColor);
            plotProbabilisticConstraints(h,voi.probCFunc,voi.dvhColor);
        else
            %Objectives
            plotNominalObjectives(h,voi.objFunc,voi.dvhColor);
            plotNominalConstraints(h,voi.cFunc,voi.dvhColor);
        end
    end
end
end

function plotNominalObjectives(h,objCell,color)
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

function plotNominalConstraints(h,constrCell,color)
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
                plot(h,dP,Vmin,'^k','MarkerFaceColor',color);
            else
                plot(h,dP,Vmax,'vk','MarkerFaceColor',color);
            end
        case 'DoseConstraints.matRad_MinMaxEUD'
            
        case 'DoseConstraints.matRad_MinMaxMeanDose'
            
        otherwise
    end
end
end

function plotProbabilisticObjectives(h,probObjCell,color)
for o = 1:numel(probObjCell)
    obj = probObjCell{o};
    objName = class(obj);
    switch objName
        case 'apm_probPieceWiseSquaredOverdose'
            plot(h,[obj.dMax obj.dMax],[0 1],'-->','Color','k','MarkerFaceColor',color);
        case 'apm_probPieceWiseSquaredUnderdose'
            plot(h,[obj.dMin obj.dMin],[0 1],'--<','Color','k','MarkerFaceColor',color);
        otherwise
    end
end


end

function plotProbabilisticConstraints(h,probConstrCell,color)
for o = 1:numel(probConstrCell)
    constr = probConstrCell{o};
    constrName = class(constr);
    
    switch constrName
        case 'apm_probDvhMinQuantileConstraint'
            dP = constr.dParam;
            Vmin = constr.minV;
            p = constr.p;
            plot(h,dP,Vmin,'^k','MarkerFaceColor',color);
        case 'apm_probDvhMaxQuantileConstraint'
            dP = constr.dParam;
            Vmax = constr.maxV;
            p = constr.p;
            plot(h,dP,Vmax,'vk','MarkerFaceColor',color);
            
        otherwise
    end
end
end