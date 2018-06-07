function apm_anatomyPlot(h,x,vois)

plotColors = apm_plotColors;

hold(h,'on');
%set(axProfile,'color','none','Layer','top','FontSize',fSize,'LineWidth',2)
%linkaxes([axHist axProfile],'y')
xlabel(h,'$x$ [mm]')
ylabel(h,'rel. dose')
%set(axProfile,'YAxisLocation','right')
%set(axProfile,'YTickLabel',[])
yLimAx1 = 1.2;
axis([x(1) x(end) 0 yLimAx1])

grid(h,'on');
box(h,'on');

%plot the vois

for v=1:numel(vois)    
    if strcmp(vois(v).type,'TARGET') || strcmp(vois(v).type,'OAR')
        
        colorGray = plotColors.anatomyBorder;
        yl = ylim(h);
        line(h,[vois(v).xL vois(v).xL],yl,'Color',colorGray);
        line(h,[vois(v).xU vois(v).xU],yl,'Color',colorGray);
        
        patch([vois(v).xL vois(v).xU vois(v).xU vois(v).xL],[yl(1) yl(1) yl(2) yl(2)],vois(v).faceColor,'Parent',h,'FaceAlpha',plotColors.anatomyAlpha);
    end
end

end

