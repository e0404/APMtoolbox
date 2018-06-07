function plotColors = apm_plotColors

allColors = colorspecs;

plotColors.nomDose = allColors.plotblue;
plotColors.expDose = allColors.plotorange;
plotColors.stdDose = allColors.plotyellow;

plotColors.target = allColors.plotblue_grayed;
plotColors.oar = allColors.plotorange_grayed;

plotColors.targetDVH = allColors.plotblue;
plotColors.oarDVH = allColors.plotorange;

plotColors.anatomyAlpha = 0.3;
plotColors.anatomyBorder = [0.5 0.5 0.5];


end

