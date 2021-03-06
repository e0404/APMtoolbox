classdef matRad_SquaredOverdosing < DoseObjectives.matRad_DoseObjective
    %MATRAD_DOSEOBJECTIVE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        name = 'Squared Overdosing';
        parameterNames = {'d^{max}'};
        parameterTypes = {'dose'};
    end
    
    properties
        parameters = {30};
        penalty = 1;
    end
    
    methods 
        %% Calculates the Objective Function value
        function fDose = computeDoseObjectiveFunction(obj,dose)                       
            % overdose : dose minus prefered dose
            overdose = dose - obj.parameters{1};
            
            % apply positive operator
            overdose(overdose<0) = 0;

            % claculate objective function
            fDose = obj.penalty/numel(dose) * (overdose'*overdose);
        end
        
        %% Calculates the Objective Function gradient
        function fDoseGrad   = computeDoseObjectiveGradient(obj,dose)
            % overdose : dose minus prefered dose
            overdose = dose - obj.parameters{1};
            
            % apply positive operator
            overdose(overdose<0) = 0;

            % calculate delta
            fDoseGrad = 2 * obj.penalty/numel(dose) * overdose;
        end
    end
    
end
