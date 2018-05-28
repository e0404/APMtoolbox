% clear environment
clear
%clf
%close all;
set(groot,'defaultlegendinterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
% set state of random engine
stream = RandStream('mt19937ar','Seed',datenum(datetime));
%stream = RandStream('mt19937ar','Seed',26689);
RandStream.setGlobalStream(stream);

allColors = colorspecs;

nomDoseColor = allColors.plotblue;
expDoseColor = allColors.plotorange;
stdDoseColor = allColors.plotyellow;

targetColor = allColors.plotblue_grayed;
oarColor = allColors.plotorange_grayed;

dvhTargetColor = allColors.plotblue;
dvhOarColor = allColors.plotorange;

fullCovariance = true;
sampleValidation = false;

% set setup error once and forever...
sigmaS = 1;
sigmaR = 2;

nFrac = 30;

latCutOff = 3;

tikzFolder = 'tikz/'; %Local tikz folder

%Spot info
spotDistance = 5;
% some parameters we need
highResFac = 5;
fSize = 14;

nVox = 100;
res = 1; %[mm]

[x,vois] = apm_createAnatomy1D(nVox,res,0.5,0.25);

xLim1 = x(1);
xLim2 = x(end);

mu = apm_setLateralBeamletGrid(x,vois,spotDistance,3*spotDistance);
n = numel(mu);

wStart     = ones(n,1);%4*rand(n,1)+2;
sigma = 2.5*ones(1,n);
mode = 'perfect';

nSamplesTotal = 100*nFrac;

samplingMethod = 'fractions'; %'fractions'; %'independent'
showWaitbar = false;
showScatterPlots = false;

Kn = @(n) sqrt((n-1)/2).*exp(gammaln((n-1)/2) - gammaln(n/2));
Vn = @(n) 2*((n-1)/2 - exp(2*gammaln(n/2) - 2*gammaln((n-1)/2)));
sampledStdAccRel = @(n) Kn(n).*sqrt(Vn(n))./sqrt(n-1);

multigauss = @(x,mu,Sigma) 1/sqrt((2*pi)^numel(x) * det(Sigma)) * exp(-0.5*(x - mu)' / Sigma * (x-mu));

if strcmp(mode,'rand')
    C_S = apm_createCovarianceMatrix(n,sigmaS,'random');
    C_R = apm_createCovarianceMatrix(n,sigmaR,'random');
elseif strcmp(mode,'perfect')
    C_S = sigmaS^2*(ones(n)+1e-10*eye(n)); % perfect correlation
    C_R = sigmaR^2*(ones(n)+1e-10*eye(n));
elseif strcmp(mode,'block')
    C1 = (ones(ceil(n/2))+1e-10*eye(ceil(n/2)));
    C2 = (ones(floor(n/2))+1e-10*eye(floor(n/2)));
    C_S  = [sigmaS*C1 zeros(ceil(n/2));zeros(floor(n/2)) sigmaS*C2];
    C_R  = [sigmaR*C1 zeros(ceil(n/2));zeros(floor(n/2)) sigmaR*C2];
elseif strcmp(mode,'uncorrelated')
    C_S = sigmaS^2*eye(n); % uncorrelated
    C_R = sigmaR^2*eye(n); 
elseif strcmp(mode,'toeplitz')
    C_S = toeplitz(logspace(0,-0.1,n)) .* ((sigmaS*ones(n,1)) * (sigmaS*ones(n,1))');
    C_R = toeplitz(logspace(0,-0.1,n)) .* ((sigmaR*ones(n,1)) * (sigmaR*ones(n,1))');
end

modeStr = '';
if ~strcmp(mode,'perfect')
    modeStr = ['_' mode];
end
    

%% Compute probabilistic dose influence

calcDose = @(dij,w) dij*w;
calcExpDose = @(edij,w) edij*w;

dose_ij = ones(nVox,1) ./ sqrt( 2*pi*ones(nVox,1)*sigma.^2) .* exp( - ( bsxfun(@minus,ones(nVox,1)*mu',x).^2 ) ./ ( 2*ones(nVox,1)*sigma.^2) );

sigmaAddSqr = sigma.^2 + diag(C_S)' + diag(C_R)';

expDose_ij = 1./sqrt( 2*pi*ones(nVox,1)*sigmaAddSqr) .* exp( - ( bsxfun(@minus,ones(nVox,1)*mu',x).^2 ) ./ ( 2*ones(nVox,1)*sigmaAddSqr) );

stdDose = NaN*x;

xStar  = mean([xLim1 xLim2])+.3*rand*(xLim2-xLim1);

varInfluenceSys = zeros(highResFac*n,n,n);
varInfluenceRand = zeros(highResFac*n,n,n);

if fullCovariance   
    covDose = NaN*ones(numel(x));

    covInfluenceSys = zeros(highResFac*n,n,highResFac*n,n);
    covInfluenceRand = zeros(highResFac*n,n,highResFac*n,n);
    for i = 1:highResFac*n
        %for l=i:highResFac*n
        for l=i:highResFac*n
            for j = 1:n
                for m = 1:n
                %for m = j:n
                    dev =  x([i l]) - mu([j m]);
                    
                    if any(abs(dev) > 5*sqrt(transpose(sigma([j m]).^2) + diag(C_S([j m],[j m])) + 1/nFrac * diag(C_R([j m],[j m]))))
                        continue;
                    end
                    %dev = mu([j k]) - x(i);
                    
                    tmp_CR = C_R([j m],[j m]);
                    
                    %Systematic Term
                    A = diag([sigma(j)^2 sigma(m)^2]) + diag(diag(tmp_CR)) + C_S([j m],[j m]);
                    fac = 1/(2*pi*sqrt(det(A)));
                    
                    vUncorr = fac * exp(-.5 * dev' / A * dev);
                    
                    %Add Random Term
                    A = A + flipud(diag(diag(flipud(tmp_CR)))); %Add off diagonals
                    fac = 1/(2*pi*sqrt(det(A)));
                    
                    vCorr = fac * exp(-.5 * dev' / A * dev);
                    
                    covInfluenceSys(i,j,l,m) = vUncorr - expDose_ij(i,j)*expDose_ij(l,m);
                    covInfluenceRand(i,j,l,m) = vCorr - vUncorr;  
                    
                    covInfluenceSys(l,m,i,j) = covInfluenceSys(i,j,l,m);
                    covInfluenceRand(l,m,i,j) = covInfluenceRand(i,j,l,m);
                end
            end
            
            %covDose(i,l) = w' * mx * w - expDose(i)*expDose(l);
            %covDose(l,i) = covDose(i,l);
        end
        %y_std(i) = sqrt(w * mx * w' - y_bar(i)^2);
        %dose_std(i) = sqrt(covDose(i,i));
        
        varInfluenceRand(i,:,:) = covInfluenceRand(i,:,i,:);
        varInfluenceSys(i,:,:) = covInfluenceSys(i,:,i,:);
        
        if (showWaitbar)
            if (i == 1 && l == 1)
                h = waitbar(i/nS_S,['Sample Treatment ' num2str(i) ' of ' num2str(nS_S)]);
            else
                waitbar(i/nS_S,h,['Sample Treatment ' num2str(i) ' of ' num2str(nS_S)]);
            end
        else
            matRad_progress(i,highResFac*n);
        end
        
    end
    if (showWaitbar)
        delete(h);
    end
    covInfluence = 1/nFrac * covInfluenceRand + covInfluenceSys;
    varInfluence = 1/nFrac * varInfluenceRand + varInfluenceSys;
    %covInfluence = tensor(covInfluence);
       
else
    c=1;
    for i = 1:highResFac*n       
        for j = 1:n
            for m = j:n
                dev =  x([i i]) - mu([j m]);
                %dev = mu([j k]) - x(i);
                
                tmp_CR = C_R([j m],[j m]);
                
                %Systematic Term
                A = diag([sigma(j)^2 sigma(m)^2]) + diag(diag(tmp_CR)) + C_S([j m],[j m]);
                fac = 1/(2*pi*sqrt(det(A)));
                
                vUncorr = fac *exp(-.5 * dev' / A * dev); 
                %mx(j,m) = (nFrac-1) * fac * exp(-.5 * dev' / A * dev);                
                %Add Random Term
                A = A + flipud(diag(diag(flipud(tmp_CR)))); %Add off diagonals
                fac = 1/(2*pi*sqrt(det(A)));
                vCorr = fac *exp(-.5 * dev' / A * dev); 
                %varInfluence(i,j,m) = (mx(j,m) + fac * exp(-.5 * dev' / A * dev )) / nFrac - expDose_ij(i,j)*expDose_ij(i,m);
                
                varInfluenceSys(i,j,m) = vUncorr - expDose_ij(i,j)*expDose_ij(i,m);
                varInfluenceSys(i,m,j) = varInfluenceSys(i,j,m);
                varInfluenceRand(i,j,m) = vCorr - vUncorr;
                varInfluenceRand(i,m,j) = varInfluenceRand(i,j,m);

            end
        end
               
        if (showWaitbar)
            if (i == 1 && l == 1)
                h = waitbar(i/nS_S,['Sample Treatment ' num2str(i) ' of ' num2str(nS_S)]);
            else
                waitbar(i/nS_S,h,['Sample Treatment ' num2str(i) ' of ' num2str(nS_S)]);
            end
        else
            matRad_progress(c,highResFac*n);
        end
        c=c+1;
    end
    if (showWaitbar)
        delete(h);
    end
       
    varInfluence = 1/nFrac * varInfluenceRand + varInfluenceSys;        
end

dij.physicalDose = dose_ij;
dij.physicalExpDose = expDose_ij;
dij.physicalCovDose = covInfluence;
dij.physicalVarDose = varInfluence;

%% Perform nominal optimization
for v=1:numel(vois)
    vois(v).objFunc = cell(0);
    vois(v).cFunc = cell(0);
    vois(v).probObjFunc = cell(0);
    vois(v).probCFunc = cell(0);
end

%obj func
%obj = 'sqDev';
obj = 'pwSqDev';
%obj = 'pwSqDevFake';
switch obj
    case 'sqDev'        
        %objFunc = @(x) (dose_ij*x - dPres)'*diag(p)*(dose_ij*x - dPres);
        %gradFunc = @(x) ((2*p.*(dose_ij*x - dPres))' * dose_ij)';
        for v=1:numel(vois)
            optFunc = DoseObjectives.matRad_SquaredDeviation;
            optFunc.parameters{1} = vois(v).dPres;
            optFunc.penalty = vois(v).p;
            vois(v).objFunc{end+1} = optFunc;        
            
            probOptFunc = apm_probLeastSquares(vois(v).p,vois(v).dPres); 
            vois(v).probObjFunc{end+1} = probOptFunc;  
        end      
    case 'pwSqDev'           
        for v=1:numel(vois)
            
            if vois(v).dObjMin == vois(v).dObjMax && vois(v).dObjMin > 0
                optFunc = DoseObjectives.matRad_SquaredDeviation;
                optFunc.parameters{1} = vois(v).dPres;
                optFunc.penalty = vois(v).p;
                vois(v).objFunc{end+1} = optFunc; 
                
                probOptFunc = apm_probLeastSquares(vois(v).p,vois(v).dPres); 
                vois(v).probObjFunc{end+1} = probOptFunc;  
            else
                if vois(v).dObjMin > 0
                    optFunc = DoseObjectives.matRad_SquaredUnderdosing;
                    optFunc.parameters{1} = vois(v).dObjMin;
                    optFunc.penalty = vois(v).p;
                    vois(v).objFunc{end+1} = optFunc;
                    
                    probOptFunc = apm_probPieceWiseSquaredUnderdose(vois(v).p,vois(v).dObjMin); 
                    vois(v).probObjFunc{end+1} = probOptFunc; 
                end
                optFunc = DoseObjectives.matRad_SquaredOverdosing;
                optFunc.parameters{1} = vois(v).dObjMax;
                optFunc.penalty = vois(v).p;
                vois(v).objFunc{end+1} = optFunc;
                
                probOptFunc = apm_probPieceWiseSquaredOverdose(vois(v).p,vois(v).dObjMax);
                vois(v).probObjFunc{end+1} = probOptFunc;
            end            
        end
    case 'pwSqDevFake'           
        for v=1:numel(vois)
            
            if vois(v).dObjMin == vois(v).dObjMax && vois(v).dObjMin > 0
                optFunc = DoseObjectives.matRad_SquaredDeviation;
                optFunc.parameters{1} = vois(v).dPres;
                optFunc.penalty = vois(v).p;
                vois(v).objFunc{end+1} = optFunc; 
                
                probOptFunc = apm_probLeastSquares(vois(v).p,vois(v).dPres); 
                vois(v).probObjFunc{end+1} = probOptFunc;  
            else
                if vois(v).dObjMin > 0
                    optFunc = DoseObjectives.matRad_SquaredUnderdosing;
                    optFunc.parameters{1} = vois(v).dObjMin;
                    optFunc.penalty = vois(v).p;
                    vois(v).objFunc{end+1} = optFunc;
                    
                    probOptFunc = apm_probFakePieceWiseSquaredUnderdose(vois(v).p,vois(v).dObjMin); 
                    vois(v).probObjFunc{end+1} = probOptFunc; 
                end
                optFunc = DoseObjectives.matRad_SquaredOverdosing;
                optFunc.parameters{1} = vois(v).dObjMax;
                optFunc.penalty = vois(v).p;
                vois(v).objFunc{end+1} = optFunc;
                
                probOptFunc = apm_probFakePieceWiseSquaredOverdose(vois(v).p,vois(v).dObjMax);
                vois(v).probObjFunc{end+1} = probOptFunc;
            end            
        end
    otherwise
        error(['Objective ''' obj ''' not implemented!']);
end
%constr = 'meanMax';
%constr = 'EUDmax';
constr = 'DVHmin';
%constr = 'DVHmax';

switch constr
    case 'DVHmin'
        dvhMinVol = 0.95;
        dvhDparam = 0.95;  
        optFunc = DoseConstraints.matRad_MinMaxDVH;
        optFunc.parameters{1} = dvhDparam;
        optFunc.parameters{2} = dvhMinVol;
        optFunc.parameters{3} = Inf;
        optFunc.referenceScalingVal = 1e-9;
        %optFunc.voxelScalingRatio = 10;
        vois(1).cFunc{end+1} = optFunc;        
        %constrFunc = @(x) dvhMinVol - vois(1).optFunc{1}.computeDoseConstraintFunction(dose_ij(vois(1).ix,:)*x);
        %constrJacob = @(x) -transpose(transpose(matRad_jacobFunc(dose_ij(vois(1).ix,:)*x,constraint,dvhDparam)) * dose_ij(vois(1).ix,:));
        %constrJacob = @(x) -transpose(transpose(vois(1).optFunc{1}.computeDoseConstraintJacobian(dose_ij(vois(1).ix,:)*x)) * dose_ij(vois(1).ix,:));
        
        dvhMinProbability = 0.15865; %1 standard deviation
        cObj = apm_probDvhMinQuantileConstraint(dvhDparam,dvhMinVol,dvhMinProbability,'normal');
        vois(1).probCFunc{end+1} = cObj;
        
    case 'EUDmin'
        eudMin = 0.99;
        %eudK = -20;
        vois(1).eudK;
        
        optFunc = DoseConstraints.matRad_MinMaxEUD;
        optFunc.parameters{1} = eudK;
        optFunc.parameters{2} = eudMin;
        optFunc.parameters{3} = Inf;
        vois(1).cFunc{end+1} = optFunc;
        
        eudMinProbability = 0.15865; %1 standard deviation
        %eudMinProbability = 0.5;
        cObj = apm_probEudMinQuantileConstraint(eudMin,eudK,eudMinProbability,'normal');
        vois(1).probCFunc{end+1} = cObj;                
    case 'minDose'
        minDose = 0.95;
        optFunc = DoseConstraints.matRad_MinMaxDose;
        optFunc.parameters{1} = minDose;
        optFunc.parameters{2} = Inf;        
        vois(1).cFunc{end+1} = optFunc;
        
        eudMinProbability = 0.05;
        cObj = apm_probEudMinQuantileConstraint(minDose,-100,eudMinProbability,'normal');
        vois(1).probCFunc{end+1} = cObj; 
    case 'DVHmax'
        dvhMaxVol = 0.3;
        dvhDparam = 0.5;
        optFunc = DoseConstraints.matRad_MinMaxDVH;
        optFunc.parameters{1} = dvhDparam;
        optFunc.parameters{2} = -Inf;
        optFunc.parameters{3} = dvhMaxVol;
        optFunc.referenceScalingVal = 1e-9;
        vois(2).cFunc{end+1} = optFunc;   
        
        dvhMaxProbability = 1-0.15865; %1 standard deviation
        cObj = apm_probDvhMaxQuantileConstraint(dvhDparam,dvhMaxVol,dvhMaxProbability,'normal');
        vois(2).probCFunc{end+1} = cObj;
        
    case 'EUDmax'
        eudMax = 0.4;
        eudK = vois(2).eudK;
        optFunc = DoseConstraints.matRad_MinMaxEUD;
        optFunc.parameters{1} = eudK;
        optFunc.parameters{2} = -Inf;
        optFunc.parameters{3} = eudMax;
        vois(2).cFunc{end+1} = optFunc;    
        
        %eudMaxProbability = 1-0.15865; %1 standard deviation
        %eudMaxProbability = 0.5; 
        eudMaxProbability = 0.95;
        cObj = apm_probEudMaxQuantileConstraint(eudMax,eudK,eudMaxProbability,'normal');
        vois(2).probCFunc{end+1} = cObj;
    case 'meanMax'
        eudMax = 0.4;
        eudK = 1;
        optFunc = DoseConstraints.matRad_MinMaxEUD;
        optFunc.parameters{1} = eudK;
        optFunc.parameters{2} = -Inf;
        optFunc.parameters{3} = eudMax;
        vois(2).cFunc{end+1} = optFunc;    
        
        %eudMaxProbability = 1-0.15865; %1 standard deviation
        %eudMaxProbability = 0.5; 
        eudMaxProbability = 0.95;
        cObj = apm_probEudMaxQuantileConstraint(eudMax,eudK,eudMaxProbability,'normal');
        vois(2).probCFunc{end+1} = cObj;
    case 'maxDose'
        maxDose = 0.75;
        optFunc = DoseConstraints.matRad_MinMaxDose;
        optFunc.parameters{1} = -Inf;
        optFunc.parameters{2} = maxDose;        
        vois(2).cFunc{end+1} = optFunc;       
        
        %maxDoseProbability = 1-0.15865;
        %maxDoseProbability = 1-0.5;
        maxDoseProbability = 0.95;
        cObj = apm_probEudMaxQuantileConstraint(maxDose,100,maxDoseProbability,'normal');
        vois(2).probCFunc{end+1} = cObj;
    otherwise       
end


%options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'Algorithm','interior-point','AlwaysHonorConstraints',true);
options = optimoptions('fmincon',...
    'Display','iter-detailed',...
    'SpecifyObjectiveGradient',true,...
    'SpecifyConstraintGradient',true,...
    'AlwaysHonorConstraints', 'bounds',...
    'ScaleProblem',true);%,...
    %'PlotFcn',{@optimplotx,@optimplotfunccount,@optimplotfval,@optimplotfirstorderopt,@optimplotconstrviolation,@optimplotstepsize,@optimplotfirstorderopt});
fungrad = @(x) apm_fminconObjFuncWrapper(x,@(x) apm_objFunc(dose_ij,x,vois),@(x) apm_objGrad(dose_ij,x,vois));
nonlcon = @(x) apm_fminconConstrFuncWrapper(x,@(x) apm_cFunc(dose_ij,x,vois),[],@(x) apm_cJacob(dose_ij,x,vois),[]);
[w,fVal] = fmincon(fungrad,wStart,[],[],[],[],zeros(n,1),Inf*ones(n,1),nonlcon,options);



%Map influence
dose = calcDose(dose_ij,w);

expDose = calcExpDose(expDose_ij,w);

covDose = apm_calcCovDose(covInfluence,w);
varDose = apm_calcVarDose(varInfluence,w);
stdDose = sqrt(varDose);

%mainFig = figure;
%set(0,'currentfigure',mainFig);
% ax 1
figure;
%axHist = axes;
%hold on
%xlabel(axHist,'rel. counts','FontSize',fSize)
%ylabel('rel. dose','FontSize',fSize)
%set(axHist,'XAxisLocation','top');
%set(gcf,'Color','w')
%set(gca,'FontSize',fSize,'LineWidth',2,'Layer','top')
%set(axHist,'ActivePositionProperty','Position')

% ax 2
axProfile = axes();%axes('Position',get(axHist,'Position'),'ActivePositionProperty','Position');
hold on
%set(axProfile,'color','none','Layer','top','FontSize',fSize,'LineWidth',2)
%linkaxes([axHist axProfile],'y')
xlabel(axProfile,'$x$ [mm]','FontSize',fSize)
ylabel('rel. dose','FontSize',fSize)
%set(axProfile,'YAxisLocation','right')
%set(axProfile,'YTickLabel',[])
yLimAx1 = 1.2;
axis([plotLim1 plotLim2 0 yLimAx1])

grid on;
box on;

%plot the vois

for v=1:nVois
    colorGray = [0.5 0.5 0.5];
    yl = ylim(axProfile);
    line(axProfile,[vois(v).xL vois(v).xL],yl,'Color',colorGray);
    line(axProfile,[vois(v).xU vois(v).xU],yl,'Color',colorGray);
    
    patch([vois(v).xL vois(v).xU vois(v).xU vois(v).xL],[yl(1) yl(1) yl(2) yl(2)],vois(v).faceColor,'Parent',axProfile,'FaceAlpha',0.3);   
end
%set(ax2,'children',flipud(get(ax2,'children')));

%figure(1)
%hold on
%plot(axProfile,x,repmat(w',nVox,1).*dose_ij,'b--')
plot(axProfile,x,dose,'Color',nomDoseColor,'LineWidth',2,'LineStyle','--')

plot(axProfile,x,expDose,'Color',expDoseColor,'LineWidth',2,'LineStyle','--')

plot(axProfile,x,stdDose,'Color',stdDoseColor,'LineWidth',2,'LineStyle','--')

grid(axProfile,'on');

%DVH
figure;
axDVH = axes;
hold(axDVH,'on');
grid(axDVH,'on');

for v = 1:nVois
    dvhs{v} = apm_DVH(dose(vois(v).ix),100,1.1);
    plot(dvhs{v}(1,:),dvhs{v}(2,:),'LineWidth',2,'LineStyle','--','Color',vois(v).dvhColor);
    %dvhsProb{v} = apm_DVH(doseProb(vois(v).ix),100,1.1);
    [expDvhs{v},stdDvhs{v}] = apm_DVHprob(expDose(vois(v).ix),covDose(vois(v).ix,vois(v).ix),100,1.1,'int_gauss');
    plot(axDVH,expDvhs{v}(1,:),expDvhs{v}(2,:),'LineWidth',2,'LineStyle','-','Color',vois(v).dvhColor);
    %axDVH.ColorOrderIndex = axDVH.ColorOrderIndex - 1;
    plot(axDVH,expDvhs{v}(1,:),expDvhs{v}(2,:)+stdDvhs{v}(2,:),'LineWidth',1,'LineStyle','-.','Color',vois(v).dvhColor);
    %axDVH.ColorOrderIndex = axDVH.ColorOrderIndex - 1;
    plot(axDVH,expDvhs{v}(1,:),expDvhs{v}(2,:)-stdDvhs{v}(2,:),'LineWidth',1,'LineStyle','-.','Color',vois(v).dvhColor);
end


switch constr
    case 'DVHmin'
        plot(dvhDparam,dvhMinVol,'^k','MarkerFaceColor','k');
    case 'DVHmax'
        plot(dvhDparam,dvhMaxVol,'vk','MarkerFaceColor','k');
    case 'minDose'
        plot([minDose minDose],[0 1],'--k<','MarkerFaceColor','k');
        [minDoseMu,minDoseStd] = apm_eudProb(expDose(vois(1).ix),covDose(vois(1).ix,vois(1).ix),-100);
        tmp_d = linspace(0,1.1,100);
        tmp_gauss = 1/sqrt(2*pi*minDoseStd^2) * exp(-0.5*(tmp_d-minDoseMu).^2 ./ minDoseStd^2);
        tmp_gauss = 0.5 * tmp_gauss / max(tmp_gauss);
        plot(tmp_d,tmp_gauss,':','Color',0.5*[1 1 1]);
        %plot([minDose minDose],[0 max(tmp_gauss)],'--k<','MarkerFaceColor','k');
        %plot([min(dose(vois(1).ix)) min(dose(vois(1).ix))],[0 max(tmp_gauss)],':','Color',0.5*[1 1 1]);
    case 'maxDose'
        plot([maxDose maxDose],[0 1],'--k<','MarkerFaceColor','k');
        [maxDoseMu,maxDoseStd] = apm_eudProb(expDose(vois(2).ix),covDose(vois(2).ix,vois(2).ix),100);
        tmp_d = linspace(0,1.1,100);
        tmp_gauss = 1/sqrt(2*pi*maxDoseStd^2) * exp(-0.5*(tmp_d-maxDoseMu).^2 ./ maxDoseStd^2);
        tmp_gauss = 0.5 * tmp_gauss / max(tmp_gauss);
        plot(tmp_d,tmp_gauss,'Color',0.5*[1 1 1]);
        %plot([maxDose maxDose],[0 max(tmp_gauss)],'--k<','MarkerFaceColor','k');
        %plot([max(dose(vois(2).ix)) max(dose(vois(2).ix))],[0 max(tmp_gauss)],':','Color',0.5*[1 1 1]);
    case 'EUDmin'
        [minEudMu,minEudStd] = apm_eudProb(expDose(vois(1).ix),covDose(vois(1).ix,vois(1).ix),eudK);
        %plot([minDose minDose],[0 1],'--k<','MarkerFaceColor','k');
        tmp_d = linspace(0,1.1,100);
        tmp_gauss = 1/sqrt(2*pi*minEudStd^2) * exp(-0.5*(tmp_d-minEudMu).^2 ./ minEudStd^2);
        tmp_gauss = 0.5 * tmp_gauss / max(tmp_gauss);
        plot(tmp_d,tmp_gauss,'Color',0.5*[1 1 1]);
        plot([eudMin eudMin],[0 max(tmp_gauss)],'-k>','MarkerFaceColor','k');
    case 'EUDmax'
        [maxEudMu,maxEudStd] = apm_eudProb(expDose(vois(2).ix),covDose(vois(2).ix,vois(2).ix),vois(2).eudK);
        
        %plot([minDose minDose],[0 1],'--k<','MarkerFaceColor','k');
        tmp_d = linspace(0,1.1,100);
        tmp_pVal = maxEudMu + maxEudStd*sqrt(2)*erfinv(2*eudMaxProbability - 1);
        [~, tmp_ix ] = min( abs( tmp_d-tmp_pVal ) );
        
        tmp_gauss = 1/sqrt(2*pi*maxEudStd^2) * exp(-0.5*(tmp_d-maxEudMu).^2 ./ maxEudStd^2);
        tmp_gauss = 0.5 * tmp_gauss / max(tmp_gauss);
        plot(tmp_d,tmp_gauss,'Color',0.5*[1 1 1]);
        plot([eudMax eudMax],[0 tmp_gauss(tmp_ix)],'-k<','MarkerFaceColor','k');  
        nomEud = apm_eud(dose(vois(2).ix),vois(2).eudK);
        plot(nomEud,0,'k*','MarkerFaceColor','k');
    case 'meanMax'
        [maxEudMu,maxEudStd] = apm_eudProb(expDose(vois(2).ix),covDose(vois(2).ix,vois(2).ix),1);
        %plot([minDose minDose],[0 1],'--k<','MarkerFaceColor','k');       
        tmp_d = linspace(0,1.1,100);
        tmp_pVal = maxEudMu + maxEudStd*sqrt(2)*erfinv(2*eudMaxProbability - 1);
        [~, tmp_ix ] = min( abs( tmp_d-tmp_pVal ) );
        tmp_gauss = 1/sqrt(2*pi*maxEudStd^2) * exp(-0.5*(tmp_d-maxEudMu).^2 ./ maxEudStd^2);
        tmp_gauss = 0.5 * tmp_gauss / max(tmp_gauss);
        plot(tmp_d,tmp_gauss,'Color',0.5*[1 1 1]);
        plot([eudMax eudMax],[0 tmp_gauss(tmp_ix)],'-k<','MarkerFaceColor','k');
        nomEud = apm_eud(dose(vois(2).ix),1);
        plot(nomEud,0,'k*','MarkerFaceColor','k');
    otherwise        
end
 
if contains(obj,'pwSqDev')
    plot([vois(2).probObjFunc{1}.dMax vois(2).probObjFunc{1}.dMax],[0 1],'--k<','MarkerFaceColor','k');
end

box(axDVH,'on');
ylim(axDVH,[0 1]);
xlim(axDVH,[0 1.1]);
xlabel(axDVH,'rel. dose');
ylabel(axDVH,'rel. volume');



%% Simulate fractions for each treatment scenario
if sampleValidation
    % 1st plot histogram around one point in backgorund
    % look at one special point and show the distribution within this graph
    % by simulating a number of treatments with fractionation
    nS_S = floor(nSamplesTotal/nFrac);
    muS    = mb_mgd(nS_S,mu,C_S)';
    %yStarR = ones(nS_R,1);
    stdMeanR = ones(nS_S,1);
    doseStarS = ones(nS_S,1);
    for s=1:nS_S
        muR = mb_mgd(nFrac,muS(:,s)',C_R)';
        tmpFracDoses = (1./ sqrt( 2*pi*ones(nFrac,1)*sigma.^2) .* exp( - ( bsxfun(@minus,xStar,muR').^2 ) ./ ( 2*ones(nFrac,1)*sigma.^2) ) )*w;
        doseStarS(s) = sum(tmpFracDoses)/nFrac;
    end
    [histX,histY]=hist(doseStarS,50);
    % normalize histogram
    binWidths = diff(histY);
    histX = (histX/binWidths(1))/nS_S;
    histAxisFac = 5;
    axis([0 max(histX)*histAxisFac 0 yLimAx1])
    barh(axHist,histY,histX,'EdgeColor','none','BarWidth',1,'FaceColor',.8*[1 1 1])
    
    
    
    nS_S = floor(nSamplesTotal/nFrac);
    
    %Sample the systematic Part
    S_S = mb_mgd(nS_S,mu,C_S)';
    doseSampleS = NaN*(ones(lowResFac*n,1)*ones(1,nS_S));
    
    %Now we obtain mean, error of mean and standard deviation for each sampled
    %treatment
    
    for i = 1:nS_S
        %Simulate nFrac fractions
        S_R = mb_mgd(nFrac,S_S(:,i)',C_R)';
        
        doseSampleF = zeros(n*lowResFac,nFrac);
        for f = 1:nFrac
            doseSampleF(:,f) = (1./ sqrt( 2*pi*ones(n*lowResFac,1)*sigma.^2) .* exp( - ( bsxfun(@minus,S_R(:,f)',xLowRes).^2 ) ./ ( 2*ones(n*lowResFac,1)*sigma.^2) ))*w;
        end
        doseSampleS(:,i) = sum(doseSampleF,2)/nFrac;
        
        %Plot the accumulated treatment dose sample
        if (showScatterPlots)
            plot(xLowRes,doseSampleS(:,i),'.','Color',[0.6 0.6 1],'MarkerSize',2);
            drawnow;
        end
        %plot(xLowRes,std(yR),'.','Color',[1 0.6 1],'MarkerSize',2)
        if (showWaitbar)
            if (i == 1)
                h = waitbar(i/nS_S,['Sample Treatment ' num2str(i) ' of ' num2str(nS_S)]);
            else
                waitbar(i/nS_S,h,['Sample Treatment ' num2str(i) ' of ' num2str(nS_S)]);
            end
        else
            matRad_progress(i,nS_S);
        end
    end
    if (showWaitbar)
        delete(h);
    end
    
    doseSampleS_mean = mean(doseSampleS,2);
    doseSampleS_std  = std(doseSampleS,[],2);
    doseSampleS_cov = cov(doseSampleS');
    
    doseSampleS_mean_stdErr = doseSampleS_std / sqrt(nS_S);
    doseSampleS_std_stdErr  = doseSampleS_std * sampledStdAccRel(nS_S);
    
    %QI
    errorbar(axProfile,xLowRes,doseSampleS_mean,doseSampleS_mean_stdErr,'rs','MarkerSize',2)
    errorbar(axProfile,xLowRes,doseSampleS_std,doseSampleS_std_stdErr,'ms','MarkerSize',2)
    %plot(xLowRes,kurtosis(yS),'y*','MarkerSize',10);
    %plot(xLowRes,skewness(yS),'k*','MarkerSize',10);
end



%Evaluate accuracy sample/analytical
if sampleValidation
    resolutionFactor = highResFac / lowResFac;
    ix = 1:numel(doseSampleS_mean);
    
    dose_mean_atSamples = expDose(ix*resolutionFactor);
    dose_std_atSamples = stdDose(ix*resolutionFactor);
    
    rmse_mean = sqrt(mean((dose_mean_atSamples - doseSampleS_mean).^2));
    rmse_std  = sqrt(mean((dose_std_atSamples - doseSampleS_std).^2));
    
    disp(['Mean RMSE: ' num2str(rmse_mean)]);
    disp(['Std RMSE: ' num2str(rmse_std)]);
end
%% Get all covariances in VOI
if fullCovariance
    ax3 = axes(figure());
    imagesc(ax3,covDose);
    colormap(hot);
    hold(ax3,'on');
    for v = 1:nVois
        vois(v).dose = dose(vois(v).ix);
        vois(v).expDose = expDose(vois(v).ix);
        vois(v).stdDose = stdDose(vois(v).ix);
        vois(v).covDose = covDose(vois(v).ix,vois(v).ix);
        ixVoiList = find(vois(v).ix);
        rectangle('Position',[ixVoiList(1),ixVoiList(2),numel(ixVoiList),numel(ixVoiList)],'LineWidth',2,'LineStyle','-','EdgeColor',[1 1 1] - 1/(v+1));
    end
    colorbar(ax3);
    %matlab2tikz([tikzFolder 'covMatrix.tikz'],'relativeDataPath','tikz/','figurehandle',ax3.Parent,'height','\figureheight','width','\figurewidth','showInfo',false,'extraAxisOptions',{'ylabel shift=-5pt','xlabel shift=-3pt'});
end

%% Optimize Probabilistic
options = optimoptions('fmincon',...
    'Display','iter-detailed',...
    'SpecifyObjectiveGradient',true,...
    'SpecifyConstraintGradient',true,...
    'AlwaysHonorConstraints', 'bounds',...
    'ScaleProblem',true);%,...
    %'PlotFcn',{@optimplotx,@optimplotfunccount,@optimplotfval,@optimplotfirstorderopt,@optimplotconstrviolation,@optimplotstepsize,@optimplotfirstorderopt});
    %@(x,optimvals,state) apm_fminconPlotGradientFcn(x,optimvals,state,probFunGrad,probNonlcon)});
    %,'OutputFcn',@apm_fminconOutputFcn);

    
probObjFunc = @(x) apm_fminconObjFuncWrapper(x,@(x) apm_probObjFunc(expDose_ij,covInfluence,x,vois),@(x) apm_probObjGrad(expDose_ij,covInfluence,x,vois));
probNonlcon = @(x) apm_fminconConstrFuncWrapper(x,@(x) apm_probCFunc(expDose_ij,covInfluence,x,vois),[],@(x) apm_probCJacob(expDose_ij,covInfluence,x,vois),[]);
[wProb,fVal] = fmincon(probObjFunc,wStart,[],[],[],[],zeros(n,1),Inf*ones(n,1),probNonlcon,options);
%[wProb,fVal] = fmincon(probFunGrad,wStart,[],[],[],[],zeros(n,1),Inf*ones(n,1),[],options);

%disp(['Final Constraint Function value: ' num2str(probConstrFunc(wProb))]);

doseProb = dose_ij*wProb;
expDoseProb = expDose_ij*wProb;
if fullCovariance
    %covDoseProb = double(ttt(tensor(covInfluence),tensor(wProb*wProb'),[4 2],[1 2])); %Product with tensor toolbox
    covDoseProb = reshape(reshape(permute(covInfluence,[1 3 4 2]),numel(x)^2,n^2) * reshape(wProb*wProb',n^2,[]),numel(x),numel(x)); %Manual product via matricication
    stdDoseProb = sqrt(diag(covDoseProb));
else
    %varDoseProb = double(ttt(tensor(varInfluence),tensor(wProb*wProb'),[3 2],[1 2])); %Product with tensor toolbox
    varDoseProb = reshape(reshape(permute(varInfluence,[1 3 2]),numel(x),n^2) * reshape(wProb*wProb',n^2,[]),numel(x),1); %Manual product via matricication
    stdDoseProb = sqrt(varDoseProb);
end

if ~isvalid(axProfile)
    figure;
    axProfile = axes(gcf);
end
plot(axProfile,x,doseProb,'Color',nomDoseColor,'LineStyle','-','LineWidth',2);
plot(axProfile,x,expDoseProb,'Color',expDoseColor,'LineStyle','-','LineWidth',2);
plot(axProfile,x,stdDoseProb,'Color',stdDoseColor,'LineStyle','-','LineWidth',2);
%plot(axProfile,x,repmat(wProb',nVox,1).*dose_ij,'b-.','LineWidth',1);

if ~exist('hAxProbDvh','var') || ~isvalid(hAxProbDvh)
    hFigProbDvh = figure;
    hAxProbDvh = axes(hFigProbDvh);
else
    clear hAxProbDvh;
end

hold(hAxProbDvh,'on');
for v = 1:nVois
    dvhsProb{v} = apm_DVH(doseProb(vois(v).ix),100,1.1);
    %[expDvhs{v},stdDvhs{v}] = apm_DVHprob(expDose(vois(v).ix),covDose(vois(v).ix,vois(v).ix),100,1.1,'int_gauss');
    plot(hAxProbDvh,dvhsProb{v}(1,:),dvhsProb{v}(2,:),'LineWidth',2,'LineStyle','--','Color',vois(v).dvhColor);
    %dvhsProb{v} = apm_DVH(doseProb(vois(v).ix),100,1.1);
    [expDvhsProb{v},stdDvhsProb{v}] = apm_DVHprob(expDoseProb(vois(v).ix),covDoseProb(vois(v).ix,vois(v).ix),100,1.1,'int_gauss');
    plot(hAxProbDvh,expDvhsProb{v}(1,:),expDvhsProb{v}(2,:),'LineWidth',2,'LineStyle','-','Color',vois(v).dvhColor);
    %hAxProbDvh.ColorOrderIndex = hAxProbDvh.ColorOrderIndex - 1;
    plot(hAxProbDvh,expDvhsProb{v}(1,:),expDvhsProb{v}(2,:)+stdDvhsProb{v}(2,:),'LineWidth',1,'LineStyle','-.','Color',vois(v).dvhColor);
    %hAxProbDvh.ColorOrderIndex = hAxProbDvh.ColorOrderIndex - 1;
    plot(hAxProbDvh,expDvhsProb{v}(1,:),expDvhsProb{v}(2,:)-stdDvhsProb{v}(2,:),'LineWidth',1,'LineStyle','-.','Color',vois(v).dvhColor);
end

box(hAxProbDvh,'on');
grid(hAxProbDvh,'on');
ylim(hAxProbDvh,[0 1]);
xlim(hAxProbDvh,[0 1.1]);
xlabel(hAxProbDvh,'rel. dose');
ylabel(hAxProbDvh,'rel. volume');

switch constr
    case 'DVHmin'
        plot(dvhDparam,dvhMinVol,'^k','MarkerFaceColor','k');
    case 'DVHmax'
        plot(dvhDparam,dvhMaxVol,'vk','MarkerFaceColor','k');
    case 'minDose'
        plot([minDose minDose],[0 1],'--k>','MarkerFaceColor','k');
        [minDoseMu,minDoseStd] = apm_eudProb(expDoseProb(vois(1).ix),covDoseProb(vois(1).ix,vois(1).ix),-100);
        tmp_d = linspace(1e-6,vois(1).dPres*1.2,100);
        tmp_gauss = 1/sqrt(2*pi*minDoseStd^2) * exp(-0.5*(tmp_d-minDoseMu).^2 ./ minDoseStd^2);
        tmp_gauss = 0.5 * tmp_gauss / max(tmp_gauss);
        plot(tmp_d,tmp_gauss,'Color',0.5*[1 1 1]);
    case 'maxDose'
        plot([maxDose maxDose],[0 max(tmp_gauss)],'--k<','MarkerFaceColor','k');
        [maxDoseMu,maxDoseStd] = apm_eudProb(expDoseProb(vois(2).ix),covDoseProb(vois(2).ix,vois(2).ix),100);
        tmp_d = linspace(0,1.1,100);
        tmp_gauss = 1/sqrt(2*pi*maxDoseStd^2) * exp(-0.5*(tmp_d-maxDoseMu).^2 ./ maxDoseStd^2);
        tmp_gauss = 0.5 * tmp_gauss / max(tmp_gauss);
        plot(tmp_d,tmp_gauss,'Color',0.5*[1 1 1]);
        %plot([maxDose maxDose],[0 max(tmp_gauss)],'--k<','MarkerFaceColor','k');
    case 'EUDmin'
        [minEudMu,minEudStd] = apm_eudProb(expDoseProb(vois(1).ix),covDoseProb(vois(1).ix,vois(1).ix),eudK);
        %plot([minDose minDose],[0 1],'--k<','MarkerFaceColor','k');
        tmp_d = linspace(0,1.1,100);
        tmp_gauss = 1/sqrt(2*pi*minEudStd^2) * exp(-0.5*(tmp_d-minEudMu).^2 ./ minEudStd^2);
        tmp_gauss = 0.5 * tmp_gauss / max(tmp_gauss);
        plot(tmp_d,tmp_gauss,'Color',0.5*[1 1 1]);
        plot([eudMin eudMin],[0 max(tmp_gauss)],'-k>','MarkerFaceColor','k');
    case 'EUDmax'
        [maxEudMu,maxEudStd] = apm_eudProb(expDoseProb(vois(2).ix),covDoseProb(vois(2).ix,vois(2).ix),eudK);
        %plot([minDose minDose],[0 1],'--k<','MarkerFaceColor','k');
        tmp_d = linspace(0,1.1,100);
        tmp_gauss = 1/sqrt(2*pi*maxEudStd^2) * exp(-0.5*(tmp_d-maxEudMu).^2 ./ maxEudStd^2);                      
        tmp_gauss = 0.5 * tmp_gauss / max(tmp_gauss);
        
        plot(tmp_d,tmp_gauss,'Color',0.5*[1 1 1]);
        
        tmp_pVal = maxEudMu + maxEudStd*sqrt(2)*erfinv(2*eudMaxProbability - 1);
        [~, tmp_ix ] = min( abs( tmp_d-tmp_pVal ) );
        
        plot([eudMax eudMax],[0 tmp_gauss(tmp_ix)],'-k<','MarkerFaceColor','k');
        
        nomEud = apm_eud(doseProb(vois(2).ix),vois(2).eudK);
        plot(nomEud,0,'k*','MarkerFaceColor','k');
        
    case 'meanMax'
        [maxEudMu,maxEudStd] = apm_eudProb(expDoseProb(vois(2).ix),covDoseProb(vois(2).ix,vois(2).ix),1);
        %plot([minDose minDose],[0 1],'--k<','MarkerFaceColor','k');
        tmp_d = linspace(0,1.1,100);
        tmp_gauss = 1/sqrt(2*pi*maxEudStd^2) * exp(-0.5*(tmp_d-maxEudMu).^2 ./ maxEudStd^2);                      
        tmp_gauss = 0.5 * tmp_gauss / max(tmp_gauss);
        
        plot(tmp_d,tmp_gauss,'Color',0.5*[1 1 1]);
        
        tmp_pVal = maxEudMu + maxEudStd*sqrt(2)*erfinv(2*eudMaxProbability - 1);
        [~, tmp_ix ] = min( abs( tmp_d-tmp_pVal ) );
        
        plot([eudMax eudMax],[0 tmp_gauss(tmp_ix)],'-k<','MarkerFaceColor','k');
        
        nomEud = apm_eud(doseProb(vois(2).ix),1);
        plot(nomEud,0,'k*','MarkerFaceColor','k');
    otherwise        
end

if contains(obj,'pwSqDev')
    plot([vois(2).probObjFunc{1}.dMax vois(2).probObjFunc{1}.dMax],[0 1],'--k<','MarkerFaceColor','k');
end

%% export
cleanfigure('handle',axProfile.Parent);
cleanfigure('handle',hAxProbDvh.Parent);
cleanfigure('handle',axDVH.Parent);
matlab2tikz([tikzFolder 'profile_' obj '_' constr '_' num2str(nFrac) 'frac' modeStr '.tikz'],'relativeDataPath','tikz/','figurehandle',axProfile.Parent,'parseStrings',false,'height','\figureheight','width','\figurewidth','showInfo',false,'extraAxisOptions',{'ylabel near ticks','xlabel near ticks'});
matlab2tikz([tikzFolder 'dvhProb_' obj '_' constr '_' num2str(nFrac) 'frac' modeStr '.tikz'],'relativeDataPath','tikz/','figurehandle',hAxProbDvh.Parent,'parseStrings',false,'height','\figureheight','width','\figurewidth','showInfo',false,'extraAxisOptions',{'ylabel near ticks','xlabel near ticks'});
matlab2tikz([tikzFolder 'dvhConv_' obj '_' constr '_' num2str(nFrac) 'frac' modeStr '.tikz'],'relativeDataPath','tikz/','figurehandle',axDVH.Parent,'parseStrings',false,'height','\figureheight','width','\figurewidth','showInfo',false,'extraAxisOptions',{'ylabel near ticks','xlabel near ticks'});
