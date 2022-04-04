% Clear environment
clear
%clf
%close all;
set(groot,'defaultlegendinterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultaxesfontsize',14);
% Set state of random engine
stream = RandStream('mt19937ar','Seed',datenum(datetime));
%stream = RandStream('mt19937ar','Seed',26689);
RandStream.setGlobalStream(stream);

plotColors = apm_plotColors;

% Anatomy
nVox = 100; %100
res = 1; %[mm]
relativeTargetSize = 0.4;
relativeOARSize = 0.2;

% Steering parameters
fullCovariance = true;
sampleValidation = true;
probabilisticOptimization = true;
separateDVHplots = false;
plotBodyDVH = false;
dvhPlotMode = 'tripleband';
%dvhPlotMode = 'beta';
probDvhMethod = 'int_gauss'; %'int_gauss', 'int_gauss_cuda', 'int_loggauss' or 'copula'
% Define copula model and marginals if 'copula' is used as probDvhMethod
copula_model = 'amh'; % 'amh' or 'gauss'
copula_marginals = cell(1,nVox);
%copula_marginals(1,1:floor(0.1*nVox)) = {'loggauss'};
%copula_marginals(1,(floor(0.1*nVox)+1):floor(0.8*nVox)) = {'shiftedbeta'};
%copula_marginals(1,(floor(0.8*nVox)+1):nVox) = {'loggauss'};
copula_marginals(:) = {'gauss'}; %'gauss', 'loggauss', 'beta', 'shiftedbeta', 'gamma', or 'gumbel'

% Irradiation geometry
spotDistance = 3; % Spot Spacing
spotWidth = 5; 
nFrac = 1; % Number of Fractions

% DVHs
nDvhBins = 500;

% Uncertainty 
sigmaS = 1;
sigmaR = 2;
correlationModel = 'perfect';

% Miscallaneous
%latCutOff = 3;
tikzFolder = 'tikz/'; %Local tikz folder


% Create Anatomy and Irradiation geometry
[x,vois] = apm_createAnatomy1D(nVox,res,relativeTargetSize,relativeOARSize);

t=tiledlayout(1,1);
axProfile = axes(t); 
%axProfile = axes('Position',get(axHist,'Position'),'ActivePositionProperty','Position');
apm_anatomyPlot(axProfile,x,vois);

spots = apm_setLateralBeamletGrid(x,vois,2.5,spotDistance,3*sqrt(sigmaR^2 + sigmaS^2));
%spots = apm_setLateralBeamletGrid(x,vois,2.5,spotDistance,3*sqrt(sigmaR^2 + sigmaS^2),0.05); %adds noise to the spot width

nSpots = numel(spots);
xLim1 = x(1);
xLim2 = x(end);

plotLim1 = x(1);
plotLim2 = x(end);

% Parameters for sample validation
lowResFac = 1; % Use one every lowResFac voxels
nLowResVox = nVox/lowResFac;
if (floor(nLowResVox) ~= nLowResVox)
    error(['ERROR: Choose a lowResFac that is an exact divider of nVox']);
end
ixLowRes = linspace(1, nLowResVox, nLowResVox) * lowResFac; % Indices of xLowRes in x
xLowRes = x(ixLowRes);
nSamplesTotal = 100;
samplingMethod = 'fractions'; %'independent'
showWaitbar = false;
showScatterPlots = false;
showAllSampDVHs = true;
showAvgSampDVHs = true;
saveSampToFile = true;

% Sample location
xStar  = mean([xLim1 xLim2])+.3*rand*(xLim2-xLim1);

% Maybe find multiple ways of setting sigmas, e.g. all equal or varying
%sigma = 2.5*ones(1,nSpots);

mu = [spots(:).mu]';
sigma = [spots(:).sigma]';
wStart = ones(nSpots,1);%4*rand(n,1)+2;

multigauss = @(x,mu,Sigma) 1/sqrt((2*pi)^numel(x) * det(Sigma)) * exp(-0.5*(x - mu)' / Sigma * (x-mu));

modeStr = '';
if ~strcmp(correlationModel,'perfect')
    modeStr = ['_' correlationModel];
end

ucm = apm_createUncertaintyModel(spots,sigmaS,sigmaR,correlationModel);
    

%% Compute probabilistic dose influence

% Calculate nominal dose influence
dose_ij = apm_calcDoseInfluenceLateral(x,spots);

% Calculate expected dose influence
expDose_ij = apm_calcExpDoseInfluenceLateral(x,spots,ucm);

% Calculate covariance influence
[covInfluenceSys,covInfluenceRand] = apm_calcCovarianceInfluenceLateral(x,spots,ucm,expDose_ij,true);

% Deduce variance influence
varInfluenceSys = zeros(nVox,nSpots,nSpots);
varInfluenceRand = zeros(nVox,nSpots,nSpots);
for i=1:nVox
    varInfluenceSys(i,:,:) = covInfluenceSys(i,:,i,:);
    varInfluenceRand(i,:,:) = covInfluenceRand(i,:,i,:);
end

%Fully fractionated (co)variance influence
covInfluence = 1/nFrac * covInfluenceRand + covInfluenceSys;
varInfluence = 1/nFrac * varInfluenceRand + varInfluenceSys;

%dij.physicalDose = dose_ij;
%dij.physicalExpDose = expDose_ij;
%dij.physicalCovDose = covInfluence;
%dij.physicalVarDose = varInfluence;


%% Nominal optimization

%Create cell arrays for objectives
for v=1:numel(vois)
    vois(v).objFunc = cell(0);
    vois(v).cFunc = cell(0);
    vois(v).probObjFunc = cell(0);
    vois(v).probCFunc = cell(0);
end

obj = 'pwSqDev';
%obj = 'sqDev';
constr = 'DVHmin';
%constr = 'DVHmax';
%constr = '';

%Set a default constraint set via helper function
vois = apm_setDefaultObjectivesAndConstraints(vois,obj,constr);

numNomConstraints = sum(cellfun(@numel,{vois(:).cFunc}));

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
[w,fVal] = fmincon(fungrad,wStart,[],[],[],[],zeros(nSpots,1),Inf*ones(nSpots,1),nonlcon,options);


%Map influence
dose = apm_calcDose(dose_ij,w);
expDose = apm_calcExpDose(expDose_ij,w);
covDose = apm_calcCovDose(covInfluence,w);
varDose = apm_calcVarDose(varInfluence,w);
stdDose = sqrt(varDose);

% Get support of shifted beta distrib. in case they are used as copula marginals
[dmin_shiftedbeta,dmax_shiftedbeta] = apm_getShiftedBetaRange(expDose,stdDose,x,vois);

% Calculate DVHs
for v = 1:numel(vois)
    disp(['Computing nominal and probabilistic DVHs for VOI ' vois(v).name '...']); 
    vois(v).nomDVH = apm_DVH(dose(vois(v).ix),nDvhBins,1.1);
    expDose(vois(v).ix)
    covDose(vois(v).ix)
    [vois(v).expDVH,vois(v).stdDVH] = apm_DVHprob(expDose(vois(v).ix),covDose(vois(v).ix,vois(v).ix),nDvhBins,1.1,probDvhMethod,copula_model,copula_marginals(vois(v).ix),dmin_shiftedbeta(vois(v).ix),dmax_shiftedbeta(vois(v).ix));
    %if (v==1)
    %    figure;
    %    ax=axes();
    %    hold(ax, 'on');
    %    plot(ax,linspace(0,1.1,nDvhBins),vois(v).nomDVH(2,:),'Color',plotColors.nomDose,'LineWidth',2);
    %    plot(ax,linspace(0,1.1,nDvhBins),vois(v).expDVH(2,:),'Color',plotColors.expDose,'LineWidth',2);
    %    plot(ax,linspace(0,1.1,nDvhBins),vois(v).stdDVH(2,:),'Color',plotColors.stdDose,'LineWidth',2);
    %    legend({'Nominal', 'Expected', 'Stdev'},'Location','southwest');
    %end
end

%% Plot result
if ~isvalid(axProfile)
    figure;
    axProfile = axes();%axes('Position',get(axHist,'Position'),'ActivePositionProperty','Position');
    apm_anatomyPlot(axProfile,x,vois);
    title('Dose Profile');
end
    
apm_profilePlot(axProfile,x,dose,expDose,stdDose,'--');

% DVH
figure;
axDVH = axes;
hold(axDVH,'on');
grid(axDVH,'on');

for v = 1:numel(vois)
    %dvhs{v} = apm_DVH(dose(vois(v).ix),100,1.1);
    %dvhsProb{v} = apm_DVH(doseProb(vois(v).ix),100,1.1);
    %[expDvhs{v},stdDvhs{v}] = apm_DVHprob(expDose(vois(v).ix),covDose(vois(v).ix,vois(v).ix),100,1.1,'int_gauss_cuda');
    if ~strcmp(vois(v).type,'BODY') || plotBodyDVH
        plot(axDVH,vois(v).nomDVH(1,:),vois(v).nomDVH(2,:),'LineWidth',0.5,'LineStyle','--','Color',vois(v).dvhColor);
        apm_plotProbDVH(axDVH,vois(v).expDVH,vois(v).stdDVH,vois(v).dvhColor,{'-','-.'},dvhPlotMode);
    end        
end

box(axDVH,'on');
ylim(axDVH,[0 1]);
xlim(axDVH,[0 1.1]);
xlabel(axDVH,'rel. dose');
ylabel(axDVH,'rel. volume');
title('Conventional optimization');

apm_plotObjConstrInDVH(axDVH,vois,false,plotBodyDVH);

% switch constr
%     case 'DVHmin'
%         plot(dvhDparam,dvhMinVol,'^k','MarkerFaceColor','k');
%     case 'DVHmax'
%         plot(dvhDparam,dvhMaxVol,'vk','MarkerFaceColor','k');
%     case 'minDose'
%         plot([minDose minDose],[0 1],'--k<','MarkerFaceColor','k');
%         [minDoseMu,minDoseStd] = apm_eudProb(expDose(vois(1).ix),covDose(vois(1).ix,vois(1).ix),-100);
%         tmp_d = linspace(0,1.1,100);
%         tmp_gauss = 1/sqrt(2*pi*minDoseStd^2) * exp(-0.5*(tmp_d-minDoseMu).^2 ./ minDoseStd^2);
%         tmp_gauss = 0.5 * tmp_gauss / max(tmp_gauss);
%         plot(tmp_d,tmp_gauss,':','Color',0.5*[1 1 1]);
%         plot([minDose minDose],[0 max(tmp_gauss)],'--k<','MarkerFaceColor','k');
%         plot([min(dose(vois(1).ix)) min(dose(vois(1).ix))],[0 max(tmp_gauss)],':','Color',0.5*[1 1 1]);
%     case 'maxDose'
%         plot([maxDose maxDose],[0 1],'--k<','MarkerFaceColor','k');
%         [maxDoseMu,maxDoseStd] = apm_eudProb(expDose(vois(2).ix),covDose(vois(2).ix,vois(2).ix),100);
%         tmp_d = linspace(0,1.1,100);
%         tmp_gauss = 1/sqrt(2*pi*maxDoseStd^2) * exp(-0.5*(tmp_d-maxDoseMu).^2 ./ maxDoseStd^2);
%         tmp_gauss = 0.5 * tmp_gauss / max(tmp_gauss);
%         plot(tmp_d,tmp_gauss,'Color',0.5*[1 1 1]);
%         plot([maxDose maxDose],[0 max(tmp_gauss)],'--k<','MarkerFaceColor','k');
%         plot([max(dose(vois(2).ix)) max(dose(vois(2).ix))],[0 max(tmp_gauss)],':','Color',0.5*[1 1 1]);
%     case 'EUDmin'
%         [minEudMu,minEudStd] = apm_eudProb(expDose(vois(1).ix),covDose(vois(1).ix,vois(1).ix),eudK);
%         plot([minDose minDose],[0 1],'--k<','MarkerFaceColor','k');
%         tmp_d = linspace(0,1.1,100);
%         tmp_gauss = 1/sqrt(2*pi*minEudStd^2) * exp(-0.5*(tmp_d-minEudMu).^2 ./ minEudStd^2);
%         tmp_gauss = 0.5 * tmp_gauss / max(tmp_gauss);
%         plot(tmp_d,tmp_gauss,'Color',0.5*[1 1 1]);
%         plot([eudMin eudMin],[0 max(tmp_gauss)],'-k>','MarkerFaceColor','k');
%     case 'EUDmax'
%         [maxEudMu,maxEudStd] = apm_eudProb(expDose(vois(2).ix),covDose(vois(2).ix,vois(2).ix),vois(2).eudK);
%         
%         plot([minDose minDose],[0 1],'--k<','MarkerFaceColor','k');
%         tmp_d = linspace(0,1.1,100);
%         tmp_pVal = maxEudMu + maxEudStd*sqrt(2)*erfinv(2*eudMaxProbability - 1);
%         [~, tmp_ix ] = min( abs( tmp_d-tmp_pVal ) );
%         
%         tmp_gauss = 1/sqrt(2*pi*maxEudStd^2) * exp(-0.5*(tmp_d-maxEudMu).^2 ./ maxEudStd^2);
%         tmp_gauss = 0.5 * tmp_gauss / max(tmp_gauss);
%         plot(tmp_d,tmp_gauss,'Color',0.5*[1 1 1]);
%         plot([eudMax eudMax],[0 tmp_gauss(tmp_ix)],'-k<','MarkerFaceColor','k');  
%         nomEud = apm_eud(dose(vois(2).ix),vois(2).eudK);
%         plot(nomEud,0,'k*','MarkerFaceColor','k');
%     case 'meanMax'
%         [maxEudMu,maxEudStd] = apm_eudProb(expDose(vois(2).ix),covDose(vois(2).ix,vois(2).ix),1);
%         plot([minDose minDose],[0 1],'--k<','MarkerFaceColor','k');       
%         tmp_d = linspace(0,1.1,100);
%         tmp_pVal = maxEudMu + maxEudStd*sqrt(2)*erfinv(2*eudMaxProbability - 1);
%         [~, tmp_ix ] = min( abs( tmp_d-tmp_pVal ) );
%         tmp_gauss = 1/sqrt(2*pi*maxEudStd^2) * exp(-0.5*(tmp_d-maxEudMu).^2 ./ maxEudStd^2);
%         tmp_gauss = 0.5 * tmp_gauss / max(tmp_gauss);
%         plot(tmp_d,tmp_gauss,'Color',0.5*[1 1 1]);
%         plot([eudMax eudMax],[0 tmp_gauss(tmp_ix)],'-k<','MarkerFaceColor','k');
%         nomEud = apm_eud(dose(vois(2).ix),1);
%         plot(nomEud,0,'k*','MarkerFaceColor','k');
%     otherwise        
% end
%  
% if contains(obj,'pwSqDev')
%     plot([vois(2).probObjFunc{1}.dMax vois(2).probObjFunc{1}.dMax],[0 1],'--k<','MarkerFaceColor','k');
% end

%% Simulate fractions for each treatment scenario
if sampleValidation
    [doseSampleS_convOpt, nS_S_convOpt] = apm_runSampleValidation(w, axDVH, expDose, stdDose, nSamplesTotal, nFrac, mu, sigma, ucm, xLowRes, ixLowRes, xStar, axProfile, t, probDvhMethod, nDvhBins, vois, showWaitbar, showScatterPlots, showAllSampDVHs, showAvgSampDVHs);
    if saveSampToFile
        filename = 'SampledDosesConvOpt.mat';
        save(filename, 'doseSampleS_convOpt', 'nS_S_convOpt', 'xLowRes', 'x', 'dose', 'expDose', 'stdDose', 'vois');
    end
end


%% Optimize Probabilistic
if (probabilisticOptimization == true)
    numProbConstraints = sum(cellfun(@numel,{vois(:).probCFunc}));
    
    % fmincon
    options = optimoptions('fmincon',...
        'algorithm','interior-point',... interior-point, sqp
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
    [wProb,fVal] = fmincon(probObjFunc,wStart,[],[],[],[],zeros(nSpots,1),Inf*ones(nSpots,1),probNonlcon,options);
    
    %minimize
    %[wProb,fVal] = minimize(probObjFunc,wStart,[],[],[],[],zeros(nSpots,1),Inf*ones(nSpots,1),probNonlcon);
    
    %ipopt
    % funcs.objective = @(x) apm_probObjFunc(expDose_ij,covInfluence,x,vois);
    % funcs.gradient = @(x) apm_probObjGrad(expDose_ij,covInfluence,x,vois);
    % funcs.constraints = @(x) apm_probCFunc(expDose_ij,covInfluence,x,vois);
    % funcs.jacobian = @(x) sparse(transpose(apm_probCJacob(expDose_ij,covInfluence,x,vois)));
    % funcs.iterfunc = @(n,f,auxdata) drawnow('update');
    % 
    % startJacob = funcs.jacobian(wStart);
    % funcs.jacobianstructure = @() sparse(ones(size(startJacob)));
    % 
    % ipoptOptions.lb = zeros(nSpots,1);
    % ipoptOptions.ub = Inf*ones(nSpots,1);
    % ipoptOptions.cu = zeros(numProbConstraints,1);
    % ipoptOptions.cl = -Inf*ones(numProbConstraints,1);
    % 
    % ipoptOptions.ipopt.hessian_approximation = 'limited-memory';
    % %options.ipopt.print_level           = 0;
    % %options.ipopt.hessian_approximation = 'limited-memory';
    % ipoptOptions.ipopt.derivative_test       = 'first-order';
    % [wProb,info] = ipopt(wStart,funcs,ipoptOptions);
    
    
    %disp(['Final Constraint Function value: ' num2str(probConstrFunc(wProb))]);
    
    
    %Map influence
    doseProb = apm_calcDose(dose_ij,wProb);
    
    expDoseProb = apm_calcExpDose(expDose_ij,wProb);
    
    covDoseProb = apm_calcCovDose(covInfluence,wProb);
    varDoseProb = apm_calcVarDose(varInfluence,wProb);
    stdDoseProb = sqrt(varDoseProb);
    
    % For shifted beta marginals
    [dmin_shiftedbeta,dmax_shiftedbeta] = apm_getShiftedBetaRange(expDoseProb,stdDoseProb,x,vois);

    % Calculate DVHs
    for v = 1:numel(vois)
        disp(['Computing nominal and probabilistic DVHs for VOI ' vois(v).name '...']); 
        vois(v).nomDVHprob = apm_DVH(doseProb(vois(v).ix),nDvhBins,1.1);
        [vois(v).expDVHprob,vois(v).stdDVHprob] = apm_DVHprob(expDoseProb(vois(v).ix),covDoseProb(vois(v).ix,vois(v).ix),nDvhBins,1.1,probDvhMethod,copula_model,copula_marginals(vois(v).ix),dmin_shiftedbeta(vois(v).ix),dmax_shiftedbeta(vois(v).ix));
    end
    
    %% Plot results
    if ~isvalid(axProfile)
        figure;
        axProfile = axes();
        %axProfile = axes('Position',get(axHist,'Position'),'ActivePositionProperty','Position');
        apm_anatomyPlot(axProfile,x,vois);
    end
    apm_profilePlot(axProfile,x,doseProb,expDoseProb,stdDoseProb,'-');
    
    if ~exist('hAxProbDvh','var') || ~isvalid(hAxProbDvh)
        hFigProbDvh = figure;
        hAxProbDvh = axes(hFigProbDvh);
    else
        clear hAxProbDvh;
    end
    
    hold(hAxProbDvh,'on');
    for v = 1:numel(vois)
        if ~strcmp(vois(v).type,'BODY') || plotBodyDVH
            plot(hAxProbDvh,vois(v).nomDVHprob(1,:),vois(v).nomDVHprob(2,:),'LineWidth',0.5,'LineStyle','--','Color',vois(v).dvhColor);
            %dvhsProb{v} = apm_DVH(doseProb(vois(v).ix),100,1.1);
        
            apm_plotProbDVH(hAxProbDvh,vois(v).expDVHprob,vois(v).stdDVHprob,vois(v).dvhColor,{'-','-.'},dvhPlotMode);
        end
        %plot(hAxProbDvh,vois(v).expDVHprob(1,:),vois(v).expDVHprob(2,:),'LineWidth',2,'LineStyle','-','Color',vois(v).dvhColor);
        %hAxProbDvh.ColorOrderIndex = hAxProbDvh.ColorOrderIndex - 1;
        %plot(hAxProbDvh,vois(v).expDVHprob(1,:),vois(v).expDVHprob(2,:)+vois(v).stdDVHprob(2,:),'LineWidth',1,'LineStyle','-.','Color',vois(v).dvhColor);
        %hAxProbDvh.ColorOrderIndex = hAxProbDvh.ColorOrderIndex - 1;
        %plot(hAxProbDvh,vois(v).expDVHprob(1,:),vois(v).expDVHprob(2,:)-vois(v).stdDVHprob(2,:),'LineWidth',1,'LineStyle','-.','Color',vois(v).dvhColor);
    end
    
    box(hAxProbDvh,'on');
    grid(hAxProbDvh,'on');
    ylim(hAxProbDvh,[0 1]);
    xlim(hAxProbDvh,[0 1.1]);
    xlabel(hAxProbDvh,'rel. dose');
    ylabel(hAxProbDvh,'rel. volume');
    title('Probabilistic optimization');
    
    apm_plotObjConstrInDVH(hAxProbDvh,vois,true,plotBodyDVH);
    
    % switch constr
    %     case 'DVHmin'
    %         plot(dvhDparam,dvhMinVol,'^k','MarkerFaceColor','k');
    %     case 'DVHmax'
    %         plot(dvhDparam,dvhMaxVol,'vk','MarkerFaceColor','k');
    %     case 'minDose'
    %         plot([minDose minDose],[0 1],'--k>','MarkerFaceColor','k');
    %         [minDoseMu,minDoseStd] = apm_eudProb(expDoseProb(vois(1).ix),covDoseProb(vois(1).ix,vois(1).ix),-100);
    %         tmp_d = linspace(1e-6,vois(1).dPres*1.2,100);
    %         tmp_gauss = 1/sqrt(2*pi*minDoseStd^2) * exp(-0.5*(tmp_d-minDoseMu).^2 ./ minDoseStd^2);
    %         tmp_gauss = 0.5 * tmp_gauss / max(tmp_gauss);
    %         plot(tmp_d,tmp_gauss,'Color',0.5*[1 1 1]);
    %     case 'maxDose'
    %         plot([maxDose maxDose],[0 max(tmp_gauss)],'--k<','MarkerFaceColor','k');
    %         [maxDoseMu,maxDoseStd] = apm_eudProb(expDoseProb(vois(2).ix),covDoseProb(vois(2).ix,vois(2).ix),100);
    %         tmp_d = linspace(0,1.1,100);
    %         tmp_gauss = 1/sqrt(2*pi*maxDoseStd^2) * exp(-0.5*(tmp_d-maxDoseMu).^2 ./ maxDoseStd^2);
    %         tmp_gauss = 0.5 * tmp_gauss / max(tmp_gauss);
    %         plot(tmp_d,tmp_gauss,'Color',0.5*[1 1 1]);
    %         %plot([maxDose maxDose],[0 max(tmp_gauss)],'--k<','MarkerFaceColor','k');
    %     case 'EUDmin'
    %         [minEudMu,minEudStd] = apm_eudProb(expDoseProb(vois(1).ix),covDoseProb(vois(1).ix,vois(1).ix),eudK);
    %         %plot([minDose minDose],[0 1],'--k<','MarkerFaceColor','k');
    %         tmp_d = linspace(0,1.1,100);
    %         tmp_gauss = 1/sqrt(2*pi*minEudStd^2) * exp(-0.5*(tmp_d-minEudMu).^2 ./ minEudStd^2);
    %         tmp_gauss = 0.5 * tmp_gauss / max(tmp_gauss);
    %         plot(tmp_d,tmp_gauss,'Color',0.5*[1 1 1]);
    %         plot([eudMin eudMin],[0 max(tmp_gauss)],'-k>','MarkerFaceColor','k');
    %     case 'EUDmax'
    %         [maxEudMu,maxEudStd] = apm_eudProb(expDoseProb(vois(2).ix),covDoseProb(vois(2).ix,vois(2).ix),eudK);
    %         %plot([minDose minDose],[0 1],'--k<','MarkerFaceColor','k');
    %         tmp_d = linspace(0,1.1,100);
    %         tmp_gauss = 1/sqrt(2*pi*maxEudStd^2) * exp(-0.5*(tmp_d-maxEudMu).^2 ./ maxEudStd^2);                      
    %         tmp_gauss = 0.5 * tmp_gauss / max(tmp_gauss);
    %         
    %         plot(tmp_d,tmp_gauss,'Color',0.5*[1 1 1]);
    %         
    %         tmp_pVal = maxEudMu + maxEudStd*sqrt(2)*erfinv(2*eudMaxProbability - 1);
    %         [~, tmp_ix ] = min( abs( tmp_d-tmp_pVal ) );
    %         
    %         plot([eudMax eudMax],[0 tmp_gauss(tmp_ix)],'-k<','MarkerFaceColor','k');
    %         
    %         nomEud = apm_eud(doseProb(vois(2).ix),vois(2).eudK);
    %         plot(nomEud,0,'k*','MarkerFaceColor','k');
    %         
    %     case 'meanMax'
    %         [maxEudMu,maxEudStd] = apm_eudProb(expDoseProb(vois(2).ix),covDoseProb(vois(2).ix,vois(2).ix),1);
    %         %plot([minDose minDose],[0 1],'--k<','MarkerFaceColor','k');
    %         tmp_d = linspace(0,1.1,100);
    %         tmp_gauss = 1/sqrt(2*pi*maxEudStd^2) * exp(-0.5*(tmp_d-maxEudMu).^2 ./ maxEudStd^2);                      
    %         tmp_gauss = 0.5 * tmp_gauss / max(tmp_gauss);
    %         
    %         plot(tmp_d,tmp_gauss,'Color',0.5*[1 1 1]);
    %         
    %         tmp_pVal = maxEudMu + maxEudStd*sqrt(2)*erfinv(2*eudMaxProbability - 1);
    %         [~, tmp_ix ] = min( abs( tmp_d-tmp_pVal ) );
    %         
    %         plot([eudMax eudMax],[0 tmp_gauss(tmp_ix)],'-k<','MarkerFaceColor','k');
    %         
    %         nomEud = apm_eud(doseProb(vois(2).ix),1);
    %         plot(nomEud,0,'k*','MarkerFaceColor','k');
    %     otherwise        
    % end
    % 
    % if contains(obj,'pwSqDev')
    %     plot([vois(2).probObjFunc{1}.dMax vois(2).probObjFunc{1}.dMax],[0 1],'--k<','MarkerFaceColor','k');
    % end
    if sampleValidation
        [doseSampleS_probOpt, nS_S_probOpt] = apm_runSampleValidation(wProb, hAxProbDvh, expDoseProb, stdDoseProb, nSamplesTotal, nFrac, mu, sigma, ucm, xLowRes, ixLowRes, xStar, axProfile, t, probDvhMethod, nDvhBins, vois, showWaitbar, showScatterPlots, showAllSampDVHs, showAvgSampDVHs);
        if saveSampToFile
            filename = 'SampledDosesProbOpt.mat';
            save(filename, 'doseSampleS_probOpt', 'nS_S_probOpt', 'xLowRes', 'x', 'doseProb', 'expDoseProb', 'stdDoseProb', 'vois');
        end
    end
end
%% Export
%cleanfigure('handle',axProfile.Parent);
%cleanfigure('handle',hAxProbDvh.Parent);
%cleanfigure('handle',axDVH.Parent);
%matlab2tikz([tikzFolder 'profile_' obj '_' constr '_' num2str(nFrac) 'frac' modeStr '.tikz'],'relativeDataPath','tikz/','figurehandle',axProfile.Parent,'parseStrings',false,'height','\figureheight','width','\figurewidth','showInfo',false,'extraAxisOptions',{'ylabel near ticks','xlabel near ticks'});
%matlab2tikz([tikzFolder 'dvhProb_' obj '_' constr '_' num2str(nFrac) 'frac' modeStr '.tikz'],'relativeDataPath','tikz/','figurehandle',hAxProbDvh.Parent,'parseStrings',false,'height','\figureheight','width','\figurewidth','showInfo',false,'extraAxisOptions',{'ylabel near ticks','xlabel near ticks'});
%matlab2tikz([tikzFolder 'dvhConv_' obj '_' constr '_' num2str(nFrac) 'frac' modeStr '.tikz'],'relativeDataPath','tikz/','figurehandle',axDVH.Parent,'parseStrings',false,'height','\figureheight','width','\figurewidth','showInfo',false,'extraAxisOptions',{'ylabel near ticks','xlabel near ticks'});
