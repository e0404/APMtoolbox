function spots = apm_setLateralBeamletGrid(x,vois,sigma,spacing,expand,relSigmaNoise)
    if nargin < 5 || isempty(expand)
        expand = spacing;
    end
    
    if nargin < 6
        relSigmaNoise = 0;
    end
    
    xL = +Inf;
    xU = -Inf;
    for v = 1:numel(vois)
        if strcmp(vois(v).type,'TARGET')
            xL = min([xL vois(v).xL]);
            xU = max([xU vois(v).xU]);
        end
    end
    
    xL = xL - expand;
    xU = xU + expand;
    
    if xL < x(1)
        xL = x(1);
    end
    
    if xU > x(end)
        xU = x(end);
    end
    
    muCell = num2cell(transpose(xL:spacing:xU));
    
    spots = struct('mu',muCell,'sigma', num2cell(sigma*(ones(numel(muCell),1) + relSigmaNoise*sigma*randn(numel(muCell),1))));
end

