function stop = apm_fminconPlotGradientFcn(x,optimValues,state,varargin)

stop = false;

fObj = varargin{1};
%{
[gradEstimate,gradEstimateErr] = gradest(fObj,x);
gradDiff = gradEstimate' - optimValues.gradient;
gradDiffRel = gradDiff./optimValues.gradient;
%}
fConst = varargin{2};
[gradEstimate,gradEstimateErr] = gradest(fConst,x);
[~,~,constrGradExact,~] = fConst(x);
gradDiff = gradEstimate' - constrGradExact;
gradDiffRel = gradDiff./constrGradExact;


switch state
    case 'iter'
        if isfield(optimValues,'gradient')
            if isscalar(optimValues.gradient)
                plotscalar(optimValues.iteration,optimValues.gradient);
            else
                plotvector(optimValues.iteration,gradDiffRel);
            end
        else
            plotvector(optimValues.iteration,optimValues.residual);
        end
end

function plotscalar(iteration,fval)
% PLOTSCALAR initializes or updates a line plot of the function value
% at each iteration.

if iteration == 0
    plotfval = plot(iteration,fval,'kd','MarkerFaceColor',[1 0 1]);
    title('Current Gradient','interp','none');
    xlabel(getString(message('MATLAB:optimfun:funfun:optimplots:LabelIteration')),'interp','none');
    set(plotfval,'Tag','optimplotfval');
    ylabel(getString(message('MATLAB:optimfun:funfun:optimplots:LabelFunctionValue')),'interp','none')
else
    plotfval = findobj(get(gca,'Children'),'Tag','optimplotfval');
    newX = [get(plotfval,'Xdata') iteration];
    newY = [get(plotfval,'Ydata') fval];
    set(plotfval,'Xdata',newX, 'Ydata',newY);
    set(get(gca,'Title'),'String',getString(message('MATLAB:optimfun:funfun:optimplots:TitleCurrentFunctionValue',sprintf('%g',fval))));
end

function plotvector(iteration,fval)
% PLOTVECTOR creates or updates a bar plot of the function values or
% residuals at the current iteration.
if iteration == 0
    xlabelText = getString(message('MATLAB:optimfun:funfun:optimplots:LabelNumberOfFunctionValues0',sprintf('%g',length(fval))));
    % display up to the first 100 values
    if numel(fval) > 100
        xlabelText = {xlabelText,getString(message('MATLAB:optimfun:funfun:optimplots:LabelShowingOnlyFirst100Values'))};
        fval = fval(1:100);
    end
    plotfval = bar(fval);
    title('Current Gradient','interp','none');
    set(plotfval,'edgecolor','none')
    set(gca,'xlim',[0,1 + length(fval)])
    xlabel(xlabelText,'interp','none')
    set(plotfval,'Tag','optimplotfval');
    ylabel(getString(message('MATLAB:optimfun:funfun:optimplots:LabelFunctionValue')),'interp','none')
else
    plotfval = findobj(get(gca,'Children'),'Tag','optimplotfval');
    % display up to the first 100 values
    if numel(fval) > 100
        fval = fval(1:100);
    end
    set(plotfval,'Ydata',fval);
end


