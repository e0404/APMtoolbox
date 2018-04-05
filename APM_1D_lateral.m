% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analytical Probabilistic Modelling of a 1D lateral dose profile
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018 Hans-Peter Wieser, Niklas Wahl, Philipp Hennig and Mark Bangert
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc,clear
addpath('utils');
SG       = @(qX,qW,qMu,qSigma)((qW/(sqrt(2*pi*(qSigma^2)))).*exp(-((qX-qMu).^2)./(2*(qSigma^2))));
MG       = @(X,MU,SIGMA)((1./((2*pi)^(3/2) .* sqrt(det(SIGMA)))).*exp(-.5 *(X-MU)*(SIGMA)^-1 *(X-MU)'));
Gauss    = @(x,mu,SqSigma) 1./(sqrt(2*pi.*SqSigma)).*exp(-((x - mu).^2./(2.*SqSigma)));
SumGauss = @(x,mu,SqSigma,w) ((1./sqrt(2*pi*ones(numel(x),1) * SqSigma') .* ...
   exp(-bsxfun(@minus,x,mu').^2 ./ (2* ones(numel(x),1) * SqSigma' ))) * w);

%% define 1D phantom geometry
Voxel.voxelSize            = 0.5;    %[mm]
targetEntry                = 50;  %[mm]
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
Voxel.presDose(Voxel.ixNT) = 0.2;
Voxel.presDose(Voxel.ixT)  = 1;
Voxel.penalty(Voxel.ixNT)  = 10;
Voxel.penalty(Voxel.ixT)   = 5000;

% define margin in [mm]
margin  = 7;

% define spot geometry and create arbritray lateral beam profile
Spot.numOfSpots = 25;
Spot.Mu         = linspace(targetEntry-margin,margin+targetExit,Spot.numOfSpots);  % define mu of each Gauss
Spot.Weights    = abs(randn([Spot.numOfSpots,1]) * 0.1 + 1.5); % define pencil beam weight
Spot.Z          = ones(Spot.numOfSpots,1);                     % set the depth dose to one
Spot.Sigma      = randn([Spot.numOfSpots 1]) * 0.1 + 2;        % define width of Gaussians normrnd(3,0.5,[NumBixel 1]);

% define uncertainties
rndSetupError  = 1;      % random error
sysSetupError  = 2;      % systematic error
numFrac        = [1 5];  % define number of fractions    note that multiple fractionation schemes can be defined e.g. [1 2 4 5]

% perfectly correlated spots
SigmaRnd   = ones(Spot.numOfSpots)  * rndSetupError^2;
SigmaSys   = ones(Spot.numOfSpots)  * sysSetupError^2;

% uncorrelated spots
% SigmaRnd   = diag(ones(Spot.numOfSpots,1))  * rndSetupError^2;
% SigmaSys   = diag(ones(Spot.numOfSpots,1))  * sysSetupError^2;

SampleRuns   = 150;  % number of samples for each individual fraction

% allocate results
dij.physicalDose    = zeros(Voxel.numOfVoxel,Spot.numOfSpots,1);
dij.physicalDoseExp = zeros(Voxel.numOfVoxel,Spot.numOfSpots,1);

for j = 1:1:Spot.numOfSpots
   dij.beamNum(j,1)            = 1;
   dij.physicalDose(:,j)    = SG(Voxel.position,Spot.Z(j),Spot.Mu(j),Spot.Sigma(j));
   dij.physicalDoseExp(:,j) = SG(Voxel.position,Spot.Z(j),Spot.Mu(j),sqrt(Spot.Sigma(j).^2 + SigmaRnd(j,j) + SigmaSys(j,j)))';
end

w0       = ones(Spot.numOfSpots,1);
options  = optimoptions('fmincon','Display','iter-detailed','GradObj','on');
func     = @(x) obFunc(dij,Voxel,x);
[w,~]    = fmincon(func,w0,[],[],[],[],zeros(Spot.numOfSpots,1),Inf*ones(Spot.numOfSpots,1),[],options);

dose     = dij.physicalDose    * w;        % nominal dose
doseExp  = dij.physicalDoseExp * w;        % expected dose

h  = figure('Color',[1 1 1]);
color = colorspecs();
%% sample scenarios from gaussion probability density
for ixFrac = 1:length(numFrac)
   
   boxTargetX = [targetEntry targetExit targetExit targetEntry ];
   boxTargetY = [0  0  1.2 1.2];
   p  = patch(boxTargetX,boxTargetY,[.8 .8 .8]);set(p,'FaceAlpha',0.35,'LineStyle','none'),hold on;box on,grid on,
   h1 = plot(Voxel.position',dose,'color',color.dkfzdB,'LineWidth',2,'DisplayName','[d]');hold on
   h2 = plot(Voxel.position',doseExp,'color',color.red	,'LineWidth',2,'DisplayName','E[d]');hold on
   
   % plot single lateral profiles
   vYsum = zeros(1,Voxel.numOfVoxel);
   for j = 1:(Spot.numOfSpots)
      vY = w(j) * SG(Voxel.position',Spot.Z(j),Spot.Mu(j),Spot.Sigma(j))';
      h3 = plot(Voxel.position',vY,'color',color.gra,'LineWidth',1,'DisplayName','lateral dose');hold on
      vYsum = vYsum + vY;
   end
   
   % determine sampling number per fraction
   numRepetition = 150;  % number of RT course samples
   mSampData     = single(zeros(Voxel.numOfVoxel,numRepetition*SampleRuns));
   cnt           = 1;
   
   for real = 1:numRepetition
      
      [U,V,S]   = eig(SigmaSys);
      deltaSys  = bsxfun(@plus,0,(U * (sqrtm(V)) * S')' *  randn(Spot.numOfSpots,SampleRuns))';
      mTotFrac  = zeros(Voxel.numOfVoxel,SampleRuns,numFrac(ixFrac));
      
      for f = 1: numFrac(ixFrac)
         
         %DeltaRnd  = randn([SampleRuns 1]) * deltaRnd;
         [U,V,S]   = eig(SigmaRnd);
         deltaRnd  = bsxfun(@plus,0,(U * (sqrtm(V)) * S')' *  randn(Spot.numOfSpots,SampleRuns))';
         vSamp     = deltaRnd + deltaSys;
         
         for ixSample = 1:SampleRuns
            
            mSample = zeros(Voxel.numOfVoxel,Spot.numOfSpots);
            for j = 1:Spot.numOfSpots
               mSample(:,j) =  SG(Voxel.position,Spot.Z(j),Spot.Mu(j)- vSamp(ixSample,j),Spot.Sigma(j))';
            end
            
            mTotFrac(:,ixSample,f) = (mSample * w);
         end
         
      end
      
      mSampData(:,cnt:SampleRuns+cnt-1) =  sum(mTotFrac,3);
      cnt = cnt + SampleRuns;
      
   end
   
   vMeanSamp = mean(mSampData,2)/numFrac(ixFrac);
   vStdSamp  = std(mSampData,1,2)/numFrac(ixFrac);
   
   plot(Voxel.position(1:5:end),vMeanSamp(1:5:end),'x','color',color.red,'LineWidth',2,'LineStyle','none'),hold on
   plot(Voxel.position(1:5:end),vStdSamp(1:5:end),'x','color',color.ora,'LineWidth',2,'LineStyle','none'),hold on
   
   % calculate variance analytically
   vStd_analy        = zeros(Voxel.numOfVoxel,1);
   mCovarianceCorr   = zeros(Voxel.numOfVoxel,Voxel.numOfVoxel,Spot.numOfSpots,Spot.numOfSpots);
   mCovarianceUnCorr = zeros(Voxel.numOfVoxel,Voxel.numOfVoxel,Spot.numOfSpots,Spot.numOfSpots);
   mOmega            = zeros(Spot.numOfSpots);
   
   for i = 1:Voxel.numOfVoxel
      
      for l = i; %i:1:Voxel.numOfVoxel
         
         mPSI_Corr_ijlm    = dij.physicalDoseExp(i,:)' * dij.physicalDoseExp(l,:);
         mPSI_UnCorr_ijlm  = dij.physicalDoseExp(i,:)' * dij.physicalDoseExp(l,:);
         
         for j = 1:Spot.numOfSpots
            for m = j:Spot.numOfSpots
               
               if SigmaRnd(j,m) > 0 || SigmaSys(j,m) > 0
                  dev_il  = [Voxel.position(i) Voxel.position(l)];
                  mu_std  = [Spot.Mu(j) Spot.Mu(m)];
                  
                  Lambda        = [Spot.Sigma(j)^2 0 ; 0 Spot.Sigma(m)^2];
                  
                  SigRndCorr    = [SigmaRnd(j,j)  SigmaRnd(j,m); SigmaRnd(m,j) SigmaRnd(m,m)];
                  SigRndUnCorr  = [SigmaRnd(j,j)  0; 0 SigmaRnd(m,m)];
                  SigSys        = [SigmaSys(j,j)  SigmaSys(j,m); SigmaSys(m,j) SigmaSys(m,m)];
                  
                  mPSI_Corr_ijlm(j,m) = (1/(2*pi*sqrt(det(Lambda + SigRndCorr + SigSys))))*exp(-0.5*(dev_il - mu_std) ...
                     * ((Lambda + SigRndCorr + SigSys)^-1) * (dev_il - mu_std)');
                  mPSI_Corr_ijlm(m,j) = mPSI_Corr_ijlm(j,m);
                  
                  mPSI_UnCorr_ijlm(j,m) = (1/(2*pi*sqrt(det(Lambda +SigRndUnCorr + SigSys))))...
                     * exp(-0.5*(dev_il - mu_std) * ((Lambda + SigRndUnCorr + SigSys)^-1) * (dev_il - mu_std)');
                  mPSI_UnCorr_ijlm(m,j) = mPSI_UnCorr_ijlm(j,m);
                  
               end
            end
         end
         
         mCovarianceCorr(i,l,:,:)   = mPSI_Corr_ijlm;
         mCovarianceUnCorr(i,l,:,:) = mPSI_UnCorr_ijlm;
         
      end
      
      PSI        = ((squeeze(mCovarianceUnCorr(i,i,:,:))*(numFrac(ixFrac)-1)) + ...
                    (squeeze(mCovarianceCorr(i,i,:,:))))/numFrac(ixFrac);
      
      PSI_Corr    = w' * squeeze(mCovarianceCorr(i,i,:,:) )  * w;
      PSI_UnCorr  = w' * squeeze(mCovarianceUnCorr(i,i,:,:)) * w;
      
      Var_Corr    = (PSI_Corr   - doseExp(i)^2);
      Var_UnCorr  = (PSI_UnCorr - doseExp(i)^2);
      VarTot      = ((Var_UnCorr*(numFrac(ixFrac)-1)) + (Var_Corr))/numFrac(ixFrac);
      
      vStd_analy(i)  =  sqrt(VarTot);
      mOmega         = mOmega + Voxel.penalty(i) * (PSI - (dij.physicalDoseExp(i,:)' * dij.physicalDoseExp(i,:)));
      
   end
   
   h4 = plot(Voxel.position,vStd_analy,'color',color.ora,'LineWidth',2);
   
   set(gcf,'units','normalized','Position',[0    0.2875    0.7500    0.700])
   xlabel('x [mm]','Interpreter','latex');
   ylabel('rel. dose','Interpreter','latex');
   title([num2str(numFrac(ixFrac)) ' fraction(s); ' num2str(sysSetupError) ' [mm] sys. error and ' num2str(rndSetupError) ' [mm] rand. error '],'Interpreter','Latex')
   
   % perform probabilistic optimization
   Voxel.presDose(Voxel.ixNT) = 0;
   w0         = 0.01 * ones(Spot.numOfSpots,1);
   options    = optimoptions('fmincon','Display','iter-detailed','GradObj','on');
   funcProb   = @(x) obFuncProb(dij,Voxel,mOmega,x);
   [wRob,~]   = fmincon(funcProb,w0,[],[],[],[],zeros(Spot.numOfSpots,1),Inf*ones(Spot.numOfSpots,1),[],options);
   
   doseRob    = dij.physicalDose    * wRob;  % nominal robust dose
   doseExpRob = dij.physicalDoseExp * wRob;  % expected robust dose
   
   % calculate robust standard deviation
   vStd_analyRob  = zeros(Voxel.numOfVoxel,1);
   
   for i = 1:Voxel.numOfVoxel
      
      PSI_Corr    = wRob' * squeeze(mCovarianceCorr(i,i,:,:))   * wRob;
      PSI_UnCorr  = wRob' * squeeze(mCovarianceUnCorr(i,i,:,:)) * wRob;
      
      Var_Corr    = (PSI_Corr   - doseExpRob(i)^2);
      Var_UnCorr  = (PSI_UnCorr - doseExpRob(i)^2);
      VarTot      = ((Var_UnCorr*(numFrac(ixFrac)-1)) + (Var_Corr))/numFrac(ixFrac);
      
      vStd_analyRob(i)  =  (sqrt(VarTot));
      
   end
   
   % plot robust std
   h5 = plot(Voxel.position,vStd_analyRob,'--','color',color.ora,'LineWidth',2);
   
end

title([num2str(numFrac) ' fraction(s); ' num2str(sysSetupError) ' [mm] sys. error and ' num2str(rndSetupError) ' [mm] rand. error '],'Interpreter','Latex')
legend([h3 h1 h2 h4 h5],{'lat. dose','[d]','E[d]','$\sigma[d]$', '$\sigma[d]$ rob.'},'Interpreter','latex');
set(gca,'FontSize',23), box on, ax = gca; ax.LineWidth = 2;



