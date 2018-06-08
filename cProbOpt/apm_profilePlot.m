function apm_profilePlot(h,x,dose,expDose,stdDose,lineStyle)

if ~isvalid(h)
    warning('No valid axis provided! Skipping...');
end

if nargin < 6
    lineStyle = '-';
end

plotColors = apm_plotColors;

hold(h,'on');

plot(h,x,dose,'Color',plotColors.nomDose,'LineWidth',2,'LineStyle',lineStyle)

plot(h,x,expDose,'Color',plotColors.expDose,'LineWidth',2,'LineStyle',lineStyle)

plot(h,x,stdDose,'Color',plotColors.stdDose,'LineWidth',2,'LineStyle',lineStyle)

grid(h,'on');
end

