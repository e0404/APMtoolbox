% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analytical Probabilistic Modelling of a 1D depth dose profile for the physical dose
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018 Hans-Peter Wieser, Niklas Wahl, Philipp Hennig and Mark Bangert
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clc,clear,%close all
addpath('utils');
fEffectToRBExD  = @(effect,alphaX,betaX)(- alphaX + sqrt(alphaX.^2 + 4 .* betaX.* effect))./(2*betaX);
SG       = @(qX,qW,qMu,qSigma)((qW/(sqrt(2*pi*(qSigma^2)))).*exp(-((qX-qMu).^2)./(2*(qSigma^2))));
MG       = @(X,MU,SIGMA)((1./((2*pi)^(3/2) .* sqrt(det(SIGMA)))).*exp(-.5 *(X-MU)*(SIGMA)^-1 *(X-MU)'));
Gauss    = @(x,mu,SqSigma) 1./(sqrt(2*pi.*SqSigma)).*exp(-((x - mu).^2./(2.*SqSigma)));
SumGauss = @(x,mu,SqSigma,w) ((1./sqrt(2*pi*ones(numel(x),1) * SqSigma') .* ...
                              exp(-bsxfun(@minus,x,mu').^2 ./ (2* ones(numel(x),1) * SqSigma' ))) * w);

load(['carbon_GenericAPM.mat']);   
%% define 1D phantom geometry
numOfBeams                 = 2;    % 1 = one beam; 2 = two opposing beams
ixTissue                   = 1;    % use alpha_x/beta_x ratio of 2
alpha_x                    = machine.data(1).alphaX(ixTissue);
beta_x                     = machine.data(1).betaX(ixTissue);
numComp                    = numel(machine.data(1).alphaDose(1).width);  
Voxel.voxelSize            = 1;    %[mm]
targetEntry                = 100;  %[mm]
targetLength               = 50;   %[mm]
targetExit                 = targetEntry + targetLength;
patientLength              = round(targetEntry + targetLength + targetEntry); %[mm]
Voxel.position             = Voxel.voxelSize/2:Voxel.voxelSize:patientLength;
Voxel.ix                   = 1:numel(Voxel.position);
Voxel.numOfVoxel           = numel(Voxel.position);
Voxel.position             = (Voxel.voxelSize/2):Voxel.voxelSize:patientLength;
Voxel.numOfVoxel           = numel(Voxel.position);
Voxel.ixNT                 = (Voxel.position < targetEntry) | (Voxel.position > targetExit);
Voxel.ixT                  = ~Voxel.ixNT ;
Voxel.presDose(Voxel.ixNT) = alpha_x * 1 + beta_x * 1.^2;
Voxel.presDose(Voxel.ixT)  = alpha_x * 3 + beta_x * 3.^2;
Voxel.penalty(Voxel.ixNT)  = 10;
Voxel.penalty(Voxel.ixT)   = 5000;

%% fill spot related variables
availablePeakPos           = [machine.data.peakPos]; 
availableranges            = [machine.data.range];
Spot.Spacing               = 1:1:numel(availablePeakPos);
availablePeakPos           = availablePeakPos(Spot.Spacing);
availableranges            = availableranges(Spot.Spacing);
margin                     = 20;

Spot.position = [];
Spot.range    = [];
Spot.ixEnergy = [];
Spot.ixBeam   = [];
radDepth      = zeros(Voxel.numOfVoxel,numOfBeams);


for i = 1:numOfBeams
   
   if i == 1
         radDepth(:,i)   = Voxel.position';
         Spot.position   = [Spot.position availablePeakPos(availablePeakPos > targetEntry-margin & availablePeakPos < targetExit+margin)];
         Spot.range      = [Spot.range availableranges(availablePeakPos     > targetEntry-margin & availablePeakPos < targetExit+margin)];
         Spot.ixEnergy   = [Spot.ixEnergy find(availablePeakPos             > targetEntry-margin & availablePeakPos < targetExit+margin)]; 
         Spot.ixBeam     = [Spot.ixBeam; i * ones(numel(availablePeakPos(availablePeakPos > targetEntry-margin & availablePeakPos < targetExit+margin)),1)];
   elseif i == 2
         radDepth(:,i)   = flip(Voxel.position');
         Spot.position   = [Spot.position availablePeakPos(availablePeakPos < targetExit+margin  & availablePeakPos > targetEntry - margin)];
         Spot.range      = [Spot.range    availableranges(availablePeakPos  < targetExit+margin  & availablePeakPos > targetEntry - margin)];
         Spot.ixEnergy   = [Spot.ixEnergy find(availablePeakPos             < targetExit+margin  & availablePeakPos > targetEntry - margin)];
         Spot.ixBeam     = [Spot.ixBeam; i * ones(numel(availablePeakPos(availablePeakPos > targetEntry-margin & availablePeakPos < targetExit+margin)),1)];
   end

end

Spot.numOfSpots            = length(Spot.position);
[Spot.numOfSpotsPerBeam,~] = histcounts(Spot.ixBeam);   

%% uncertainty definition
sampleRuns                 = 5000;
sysRangeError              = 0.035;    % 3.5% relative systematic range error
rndRangeError              = 1;        % 1mm absolute random range error
sysBioError                = 0.15;     % 15% error in ion response LQM parameter

% define covariance matrix skelet 
if numOfBeams == 1
   mCovDepthSkelet = blkdiag(ones(Spot.numOfSpotsPerBeam(1)));
else
   mCovDepthSkelet = blkdiag(ones(Spot.numOfSpotsPerBeam(1)),ones(Spot.numOfSpotsPerBeam(2)));
end

% obtain mean and covariance of biological uncertainties
[vMeanAlphaDose,mCovAlphaDose,vMeanSqrtBetaDose,mCovSqrtBetaDose,mCovBioTot] = getBioCovariance(machine,Spot,sysBioError,ixTissue,numComp);

%% calculate nominal and expected dose influence data
dij.physicalDose     = spalloc(Voxel.numOfVoxel,Spot.numOfSpots,1);
dij.physicalDoseExp  = spalloc(Voxel.numOfVoxel,Spot.numOfSpots,1);
dij.mAlphaDose       = spalloc(Voxel.numOfVoxel,Spot.numOfSpots,1);
dij.mAlphaDoseExp    = spalloc(Voxel.numOfVoxel,Spot.numOfSpots,1);
dij.mSqrtBetaDose    = spalloc(Voxel.numOfVoxel,Spot.numOfSpots,1);
dij.mSqrtBetaDoseExp = spalloc(Voxel.numOfVoxel,Spot.numOfSpots,1);
vSpotRange           = zeros(Spot.numOfSpots,1); 

for j = 1:Spot.numOfSpots
   
   dij.beamNum(j,1) = Spot.ixBeam(j);
   baseEntry        = machine.data(Spot.ixEnergy(j));
   deltaZ           = rndRangeError^2;
   deltaZ           = deltaZ + baseEntry.range^2 * sysRangeError^2;
   vSpotRange(j)    = baseEntry.range;
     
   dij.physicalDose(:,j)     =  SumGauss(radDepth(:,Spot.ixBeam(j)),baseEntry.Z.mean,(baseEntry.Z.width).^2,          baseEntry.Z.weight);
   dij.physicalDoseExp(:,j)  =  SumGauss(radDepth(:,Spot.ixBeam(j)),baseEntry.Z.mean,(baseEntry.Z.width).^2 + deltaZ, baseEntry.Z.weight);  
   
   dij.mAlphaDose(:,j)       =  SumGauss(radDepth(:,Spot.ixBeam(j)),baseEntry.alphaDose(ixTissue).mean,(baseEntry.alphaDose(ixTissue).width).^2,          baseEntry.alphaDose(ixTissue).weight);
   dij.mAlphaDoseExp(:,j)    =  SumGauss(radDepth(:,Spot.ixBeam(j)),baseEntry.alphaDose(ixTissue).mean,(baseEntry.alphaDose(ixTissue).width).^2 + deltaZ, baseEntry.alphaDose(ixTissue).weight);  
 
   dij.mSqrtBetaDose(:,j)    =  SumGauss(radDepth(:,Spot.ixBeam(j)),baseEntry.SqrtBetaDose(ixTissue).mean,(baseEntry.SqrtBetaDose(ixTissue).width).^2,          baseEntry.SqrtBetaDose(ixTissue).weight);
   dij.mSqrtBetaDoseExp(:,j) =  SumGauss(radDepth(:,Spot.ixBeam(j)),baseEntry.SqrtBetaDose(ixTissue).mean,(baseEntry.SqrtBetaDose(ixTissue).width).^2 + deltaZ, baseEntry.SqrtBetaDose(ixTissue).weight);  
    
end
    
w0       = 0.2*ones(Spot.numOfSpots,1);
options  = optimoptions('fmincon','Display','iter-detailed','GradObj','on');
options.OptimalityTolerance = 1e-8;
options.MaxIterations       = 3000;
func     = @(x) obFunc(dij,Voxel,x);
[w,~]    = fmincon(func,w0,[],[],[],[],zeros(Spot.numOfSpots,1),Inf*ones(Spot.numOfSpots,1),[],options);

effect    = (dij.mAlphaDose    * w)  + (dij.mSqrtBetaDose * w).^2;        % nominal dose
effectExp = (dij.mAlphaDoseExp * w)  + (dij.mSqrtBetaDose * w).^2;        % expected dose

%% calculate standard deviation using APM
std_d            = zeros(Voxel.numOfVoxel,1);
mCovariance      = zeros(Voxel.numOfVoxel,Voxel.numOfVoxel,Spot.numOfSpots,Spot.numOfSpots);
mCovariance_phys = zeros(Voxel.numOfVoxel,Voxel.numOfVoxel,Spot.numOfSpots,Spot.numOfSpots);
numComp          = length(machine.data(1).Z.weight);
 
mSysCovRadDept   = sysRangeError^2  * (vSpotRange * vSpotRange') .* mCovDepthSkelet;
mRndCovRadDepth  = rndRangeError^2 .*  mCovDepthSkelet;
mSysCovBio       = ones(Spot.numOfSpots);

mOmega = zeros(Spot.numOfSpots);

f = waitbar(0,'Please wait...');
for i = 1:1:Voxel.numOfVoxel
   
   for l = i %i:1:Voxel.numOfVoxel
      
      PSI_ijlm      = dij.mAlphaDoseExp(i,:)' * dij.mAlphaDoseExp(l,:);
      PSI_ijlm_phys = dij.mAlphaDoseExp(i,:)' * dij.mAlphaDoseExp(l,:);
      
      mW_Cov_jo   = zeros(numComp); 
      
      for j = 1:Spot.numOfSpots
         
         baseEntryJ = machine.data(Spot.ixEnergy(j));
         SigmaSq_J  = baseEntryJ.alphaDose(ixTissue).width.^2;
         Weight_J   = baseEntryJ.alphaDose(ixTissue).weight;
         Dev_J      = radDepth(i,Spot.ixBeam(j)) - baseEntryJ.alphaDose(ixTissue).mean;
         
         for m = j:Spot.numOfSpots
            if mSysCovRadDept(j,m) > 0 || mRndCovRadDepth(j,m) > 0 || mSysCovBio(j,m) > 0  % check for correlation
               
               baseEntryM = machine.data(Spot.ixEnergy(m));
               SigmaSq_M  = baseEntryM.alphaDose(ixTissue).width.^2;
               Weight_M   = baseEntryM.alphaDose(ixTissue).weight;
               Dev_M      = radDepth(l,Spot.ixBeam(m)) - baseEntryM.alphaDose(ixTissue).mean;
               
               vLaSi11 = SigmaSq_J + mSysCovRadDept(j,j) + mRndCovRadDepth(j,j);
               vLaSi22 = SigmaSq_M + mSysCovRadDept(m,m) + mRndCovRadDepth(m,m);
               vLaSi12 = mSysCovRadDept(j,m);
               vLaSi21 = mSysCovRadDept(m,j);
               
               mW_CovAlphaDose_jm = zeros(numComp);
               PSI_ijlm_phys(j,m) = calcSecRangeMom(vLaSi11,vLaSi22,vLaSi12,vLaSi21,Dev_J,Dev_M,Weight_J,Weight_M,mW_CovAlphaDose_jm);
               PSI_ijlm_phys(m,j) = PSI_ijlm_phys(j,m);  
               
               if mSysCovBio(j,m) > 0
                  mW_CovAlphaDose_jm = mCovAlphaDose(j*numComp-numComp+1:j*numComp,m*numComp-numComp+1:m*numComp);
               end
               
               PSI_ijlm(j,m) = calcSecRangeMom(vLaSi11,vLaSi22,vLaSi12,vLaSi21,Dev_J,Dev_M,Weight_J,Weight_M,mW_CovAlphaDose_jm);
               PSI_ijlm(m,j) = PSI_ijlm(j,m);  
               
            end
         end
      end
      mCovariance(i,l,:,:)      = PSI_ijlm;
      mCovariance(l,i,:,:)      = PSI_ijlm;
      mCovariance_phys(i,l,:,:) = PSI_ijlm_phys;
      mCovariance_phys(l,i,:,:) = PSI_ijlm_phys;
   end
   
   mOmega = mOmega + Voxel.penalty(i) * (squeeze(mCovariance_phys(i,i,:,:)) - dij.mAlphaDoseExp(i,:)' * dij.mAlphaDoseExp(i,:));
   
   std_d(i) = sqrt((w'* squeeze(mCovariance(i,i,:,:)) * w) - (dij.mAlphaDoseExp(i,:) * w)^2);
   waitbar(i/Voxel.numOfVoxel,f,'Calculating (co)-variance...');
   
end
close(f)

% perform probabilistic optimization
Voxel.presDose(Voxel.ixNT) = 0;
w0         = 0.01 * ones(Spot.numOfSpots,1);
funcProb   = @(x) obFuncProb(dij,Voxel,mOmega,x);
[wRob,~]   = fmincon(funcProb,w0,[],[],[],[],zeros(Spot.numOfSpots,1),Inf*ones(Spot.numOfSpots,1),[],options);

effectRob    = (dij.mAlphaDose    * wRob)  + (dij.mSqrtBetaDose * wRob).^2;        % nominal dose
effectExpRob = (dij.mAlphaDoseExp * wRob)  + (dij.mSqrtBetaDose * wRob).^2;        % expected dose

%% calculate robust std
std_d_rob = zeros(Voxel.numOfVoxel,1);

f = waitbar(0,'Please wait...');
for i = 1:1:Voxel.numOfVoxel
   
   std_d_rob(i) = (sqrt((wRob'* squeeze(mCovariance(i,i,:,:)) * wRob) - (dij.mAlphaDoseExp(i,:) * wRob)^2));
   waitbar(i/Voxel.numOfVoxel,f,'Calculating (co)-variance...');
   
end
close(f)

%% sampling non robust and robust scenarios 
[U,V,S]       = eig(mSysCovRadDept);
vSampRangeSys = bsxfun(@plus,0,(U *real(sqrtm(V))*S')' * randn(numel(vSpotRange),sampleRuns,1))';

[U,V,S]       = eig(mRndCovRadDepth);
vSampRangeRnd = bsxfun(@plus,0,(U *real(sqrtm(V))*S')' * randn(numel(vSpotRange),sampleRuns,1))';

[U,V,S]    = eig(mCovBioTot); 
mSampBioWeigh = bsxfun(@plus,[vMeanAlphaDose; vMeanSqrtBetaDose],(U * real(sqrtm(V)) * S')' *  randn(Spot.numOfSpots*numComp*2,sampleRuns))';
  
h = waitbar(0,'Please wait...');
mSampDose    = zeros(Voxel.numOfVoxel,sampleRuns);
mSampDoseRob = zeros(Voxel.numOfVoxel,sampleRuns);
Offset       = Spot.numOfSpots*numComp;

for ixSample = 1:sampleRuns
   
    mSampleAlphaDose = zeros(Voxel.numOfVoxel,Spot.numOfSpots);
    mSampleSqrtBetaDose = zeros(Voxel.numOfVoxel,Spot.numOfSpots);
    
    for j = 1:Spot.numOfSpots
        baseEntry = machine.data(Spot.ixEnergy(j));                                     
        mSampleAlphaDose(:,j)    =   SumGauss(radDepth(:,Spot.ixBeam(j)),baseEntry.alphaDose(ixTissue).mean + vSampRangeSys(ixSample,j)    + vSampRangeRnd(ixSample,j) ,...
                                              (baseEntry.alphaDose(ixTissue).width).^2 ,mSampBioWeigh(ixSample,j*numComp-numComp+1:j*numComp)'); 
        mSampleSqrtBetaDose(:,j) =   SumGauss(radDepth(:,Spot.ixBeam(j)),baseEntry.SqrtBetaDose(ixTissue).mean + vSampRangeSys(ixSample,j) + vSampRangeRnd(ixSample,j) ,...
                                              (baseEntry.SqrtBetaDose(ixTissue).width).^2 ,mSampBioWeigh(ixSample,j*numComp-numComp+1+Offset:(j*numComp)+Offset)'); 
    end

    effectSamp    = (mSampleAlphaDose  * w   )  + (mSampleSqrtBetaDose * w   ).^2; 
    effectSampRob = (mSampleAlphaDose  * wRob)  + (mSampleSqrtBetaDose * wRob).^2; 
   
    mSampDose(:,ixSample)    = fEffectToRBExD(effectSamp,alpha_x,beta_x);
    mSampDoseRob(:,ixSample) = fEffectToRBExD(effectSampRob,alpha_x,beta_x); 
    waitbar(ixSample/sampleRuns,h,'Sampling ...');
    
end

close(h)
RBExD_Exp_samp    = mean(mSampDose,2);
RBExD_Std_samp    = std(mSampDose,1,2);
RBExD_ExpRob_samp = mean(mSampDoseRob,2);
RBExD_StdRob_samp = std(mSampDoseRob,1,2);

%% error propagation
RBExD     = fEffectToRBExD(effect,alpha_x,beta_x);
RBExD_Rob = fEffectToRBExD(effect,alpha_x,beta_x);

[RBExD_Exp,RBExD_Std]       = MomEffectToMomRBExDose(effectExp,    std_d ,    alpha_x * ones(Voxel.numOfVoxel,1),beta_x * ones(Voxel.numOfVoxel,1));
[RBExD_ExpRob,RBExD_StdRob] = MomEffectToMomRBExDose(effectExpRob, std_d_rob ,alpha_x * ones(Voxel.numOfVoxel,1),beta_x * ones(Voxel.numOfVoxel,1));

%% plot everything
color = colorspecs();
boxTargetX = [targetEntry  targetEntry+targetLength  targetEntry+targetLength       targetEntry ];
boxTargetY = [0            0                         4 4];
sampSteps  = round(Voxel.voxelSize^-1)*2;

figure,set(gcf,'Color',[1 1 1],'Units','normalized'), hold on
subplot(121),p  = patch(boxTargetX,boxTargetY,[.8 .8 .8]);set(p,'FaceAlpha',0.35,'LineStyle','none'),hold on;box on,
for j = 1:Spot.numOfSpots    
    baseEntryJ = machine.data(Spot.ixEnergy(j));
    subplot(121),plot(Voxel.position, SumGauss(radDepth(:,Spot.ixBeam(j)),baseEntryJ.Z.mean,...
         (baseEntryJ.Z.width).^2,baseEntryJ.Z.weight)*w(j),'color',color.gra); hold on
end

subplot(121);hNom     = plot(Voxel.position,RBExD,'k','LineWidth',2);hold on 
subplot(121);hExp     = plot(Voxel.position,RBExD_Exp,'color',color.dre,'LineWidth',2); hold on
subplot(121);hExpSamp = plot(Voxel.position(1:sampSteps:end),RBExD_Exp_samp(1:sampSteps:end),'x','color',color.dre,'LineWidth',2); hold on
subplot(121);hStd     = plot(Voxel.position,RBExD_Std,'color',color.ora,'LineWidth',2); hold on
subplot(121);hStdSamp = plot(Voxel.position(1:sampSteps:end),RBExD_Std_samp(1:sampSteps:end),'x','color',color.ora,'LineWidth',2); hold on

set(gca,'YLim',[0 4]),set(gca,'XLim',[0 Voxel.position(end)])
xlabel('[mm]','Interpreter','Latex'),ylabel('dose [Gy(RBE)]','Interpreter','Latex'),
grid on, set(gca,'FontSize',16),
title('SOBP C12 - conventional optimization','Interpreter','Latex','FontSize',20);
legend([hNom hExp hStd],{'[RBExD]','E[RBExD] ','$\sigma[RBExD]$  '},'Interpreter','Latex','FontSize',20);
ax = gca; ax.LineWidth = 2;

subplot(122),p  = patch(boxTargetX,boxTargetY,[.8 .8 .8]);set(p,'FaceAlpha',0.35,'LineStyle','none'),hold on;box on,
for j = 1:Spot.numOfSpots    
    baseEntryJ = machine.data(Spot.ixEnergy(j));
    subplot(122),plot(Voxel.position, SumGauss(radDepth(:,Spot.ixBeam(j)),baseEntryJ.Z.mean,...
         (baseEntryJ.Z.width).^2,baseEntryJ.Z.weight)*wRob(j),'color',color.gra); hold on
end

subplot(122);hNom     = plot(Voxel.position,RBExD_Rob,'k','LineWidth',2);hold on 
subplot(122);hExp     = plot(Voxel.position,RBExD_ExpRob,'color',color.dre,'LineWidth',2); hold on
subplot(122);hExpSamp = plot(Voxel.position(1:sampSteps:end),RBExD_ExpRob_samp(1:sampSteps:end),'x','color',color.dre,'LineWidth',2); hold on
subplot(122);hStd     = plot(Voxel.position,RBExD_StdRob,'color',color.ora,'LineWidth',2); hold on
subplot(122);hStdSamp = plot(Voxel.position(1:sampSteps:end),RBExD_StdRob_samp(1:sampSteps:end),'x','color',color.ora,'LineWidth',2); hold on

set(gca,'YLim',[0 4]),set(gca,'XLim',[0 Voxel.position(end)])
xlabel('[mm]','Interpreter','Latex'),ylabel('dose [Gy(RBE)]','Interpreter','Latex'),
grid on,set(gca,'FontSize',16),
title('SOBP C12 - probabilistic optimization','Interpreter','Latex','FontSize',20);
legend([hNom hExp hStd],{'[RBExD]','E[RBExD]','$\sigma[RBExD]$  '},'Interpreter','Latex','FontSize',20);
ax = gca; ax.LineWidth = 2;


