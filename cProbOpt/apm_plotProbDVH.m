function apm_plotProbDVH(h,expDvh,stdDvh,color,styles,mode)

if ~isvalid(h)
    warning('No valid axis provided! Skipping...');
end

if nargin < 6
    mode = 'none';
end

d = expDvh(1,:);
mu = expDvh(2,:);
sig = stdDvh(2,:);

plot(h,d,mu,'Color',color,'LineStyle',styles{1},'LineWidth',2);
%plot(h,d,mu + sig,'Color',color,'LineStyle',styles{2},'LineWidth',1);
%plot(h,d,mu - sig,'Color',color,'LineStyle',styles{2},'LineWidth',1);

switch mode
    case 'band'
        fill(h,[d';flipud(d')],[(mu-sig)';flipud((mu+sig)')],color,'FaceAlpha',0.35,'LineStyle','none');
    case 'doubleband'
        fill(h,[d';flipud(d')],[(mu-sig)';flipud((mu+sig)')],color,'FaceAlpha',0.25,'LineStyle','none');
        fill(h,[d';flipud(d')],[(mu-2*sig)';flipud((mu+2*sig)')],color,'FaceAlpha',0.2,'LineStyle','none');        
    case 'tripleband'
        fill(h,[d';flipud(d')],[(mu-sig)';flipud((mu+sig)')],color,'FaceAlpha',0.2,'LineStyle','none');
        fill(h,[d';flipud(d')],[(mu-2*sig)';flipud((mu+2*sig)')],color,'FaceAlpha',0.15,'LineStyle','none');
        fill(h,[d';flipud(d')],[(mu-3*sig)';flipud((mu+3*sig)')],color,'FaceAlpha',0.1,'LineStyle','none');
        
    case 'normal'
        dvSpace = linspace(0,1,100);
        dvMat = repmat(dvSpace',1,numel(d));
        muMat = repmat(mu,100,1);
        sigMat = repmat(sig,100,1);
        p = arrayfun(@normpdf,dvMat,muMat,sigMat);
        pMax = max(p);
        pMax = repmat(pMax,100,1);
        alphaMask = p./pMax * 0.5;
        
        cData = ones([size(p) 3]);
        cData(:,:,1) = color(1);
        cData(:,:,2) = color(2);
        cData(:,:,3) = color(3);
        
        hIm = imagesc('XData',d,'YData',dvSpace,'CData',cData,'AlphaData',alphaMask,'Parent',h);
        %uistack(hIm,'bottom');
        
    case 'beta'
        dvSpace = linspace(0,1,100);
        %dvSpace = dvSpace(2:end-1);
        dvMat = repmat(dvSpace',1,numel(d));
        
        
        %parametric transform
        [a,b] = apm_transformMeanVarianceToBetaParameters(mu,sig);
        
        aMat = repmat(a,numel(dvSpace),1);
        bMat = repmat(b,numel(dvSpace),1);
        
        
        
        p = arrayfun(@betapdf,dvMat,aMat,bMat);
        p(isinf(p)) = 1; %Hack for infinity at 1
        
        pMax = max(p);
        pMax = repmat(pMax,numel(dvSpace),1);
        alphaMask = p./pMax * 0.5;
        
        cData = ones([size(p) 3]);
        cData(:,:,1) = color(1);
        cData(:,:,2) = color(2);
        cData(:,:,3) = color(3);
        
        imagesc('XData',d,'YData',dvSpace,'CData',cData,'AlphaData',alphaMask,'Parent',h);
        
    case 'truncated'
        dvSpace = linspace(0,1,100);
        dvMat = repmat(dvSpace',1,numel(d));
        muMat = repmat(mu,100,1);
        sigMat = repmat(sig,100,1);
        p = arrayfun(@truncatedNormpdf,dvMat,muMat,sigMat);
        pMax = max(p);
        pMax = repmat(pMax,100,1);
        alphaMask = p./pMax * 0.5;
        
        cData = ones([size(p) 3]);
        cData(:,:,1) = color(1);
        cData(:,:,2) = color(2);
        cData(:,:,3) = color(3);
        
        hIm = imagesc('XData',d,'YData',dvSpace,'CData',cData,'AlphaData',alphaMask,'Parent',h);
        
    otherwise
        p = str2double(mode);
        if p < 1 && p > 0
            q = norminv(p,expDvh,stdDvh);
            fill(h,[d';flipud(d')],[(mu-sig)';flipud((mu+sig)')],color,'FaceAlpha',0.2,'LineStyle','none');
            fill(h,[d';flipud(d')],[(mu-2*sig)';flipud((mu+2*sig)')],color,'FaceAlpha',0.15,'LineStyle','none');
            fill(h,[d';flipud(d')],[(mu-3*sig)';flipud((mu+3*sig)')],color,'FaceAlpha',0.1,'LineStyle','none');
            plot(h,d,q,'Color',color,'LineStyle',styles{2},'LineWidth',2);
        else
            
        end                
end

end

function p = truncatedNormpdf(x,mu,sig)
    pd = makedist('Normal','mu',mu,'sigma',sig);
    
    pd = truncate(pd,0,1);
    p = pd.pdf(x);
end
