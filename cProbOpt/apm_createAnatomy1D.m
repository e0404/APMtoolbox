function [x,vois] = apm_createAnatomy1D(nVox,res,relTargetSize,relOarSize)

width = nVox * res;
x = linspace(-width/2,width/2,nVox)';

absTargetSize = relTargetSize * abs(width);
absOarSize = relOarSize * abs(width);

allColors = colorspecs;
targetColor = allColors.plotblue_grayed;
oarColor = allColors.plotorange_grayed;
dvhTargetColor = allColors.plotblue;
dvhOarColor = allColors.plotorange;


vois = struct();
%Target
vois(1).xL = -absTargetSize/2;
vois(1).xU = +absTargetSize/2;
vois(1).dPres = 1;
vois(1).dObjMin = 1;
vois(1).dObjMax = 1;
vois(1).p = 5000;
vois(1).name = 'TARGET';
vois(1).type = 'TARGET';
vois(1).faceColor = targetColor;
vois(1).dvhColor = dvhTargetColor;
vois(1).eudK = -20;

%OAR
vois(2).xL = absTargetSize/2;
vois(2).xU = absTargetSize/2 + absOarSize;
vois(2).dPres = 0;
vois(2).dObjMin = 0;
vois(2).dObjMax = 0.6;%0.4; %0.6; %0.5;
vois(2).p = 200;
vois(2).name = 'OAR';
vois(2).type = 'OAR';
vois(2).faceColor = oarColor;
vois(2).dvhColor = dvhOarColor;
vois(2).eudK = 5.5;

nVois = numel(vois);
allVoiIx = [];

for v = 1:nVois
    %vois(v).ixLowRes = and(xLowRes <= vois(v).xU,xLowRes >= vois(v).xL);
    %vois(v).xLowRes = x(vois(v).ixLowRes);
    vois(v).ix = and(x <= vois(v).xU,x >= vois(v).xL);
    allVoiIx = [allVoiIx vois(v).ix];
    vois(v).x = x(vois(v).ix);
      
    vois(v).objFunc = cell(0,0);
    vois(v).cFunc = cell(0,0);
end

allVoiIx = any(allVoiIx,2);
%allVoiX = x(allVoiIx);
%bodyIx = 1:numel(x);

vois(3).dPres = 0;
vois(3).dObjMin = 0;
vois(3).dObjMax = 0.6;
vois(3).p = 1;
vois(3).name = 'BODY';
vois(3).type = 'BODY';
vois(3).ix = ~allVoiIx;
vois(3).x = x(vois(3).ix);
vois(3).objFunc = cell(0,0);
vois(3).cFunc = cell(0,0);
vois(3).eudK = 1;
vois(3).faceColor = 'none';
vois(3).dvhColor = [0 0 0];

end

