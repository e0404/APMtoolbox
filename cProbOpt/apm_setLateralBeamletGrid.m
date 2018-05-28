function muPos = apm_setLateralBeamletGrid(x,vois,spacing,expand)
    if nargin < 4
        expand = spacing;
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
    
    muPos = transpose(xL:spacing:xU);    
end

