 clc,clear, close all
 addpath([pwd filesep 'utils']);
 
 %% load the dataset that should be fitted
load('protons_Generic.mat')                                        % load matRad machine file
vNumComp          = [10];                                          % number of components that should be fitted - you can use multiple number of components to see what works best
vEnergyIx         = [1:1:length(machine.data)];                    % determine specific energu indices for fitting
Param.Optimizer   = {'minimize','fmincon'};                        % 'IPOPT' 'minimize','fmincon','fminsearch'; you can also combine multiple subsequent optimizations
                                                                   % if IPOPT is used make sure the 'optimization' folder of matRad is accessable and added to the path 
Param.Types       = {'Z'};                                         % determine which component should be fitted 'SqrtBetaDose';%'alphaDose';% 'Z' {'Z','alphaDose','SqrtBetaDose'};
Param.radMod      = machine.meta.radiationMode; 
Param.vEnergies   = [machine.data.energy];   
Param.bioIx       = [1];                                           % determines reference tissue for alpha and beta curves (carbon ion) machine.data(ix).alphaDose(:,param.bioIx)
                                                                   % set Param.bioIx = 1 if physical ddd should is fitted
visBool           = 1;                                             % plot fits on each iteration 
visAllBool        = 1;                                             % plot fits after successfull fitting

% function handle for calculating depth doses
SumGauss = @(x,mu,SqSigma,w) ((1./sqrt(2*pi*ones(numel(x),1) * SqSigma') .*  exp(-bsxfun(@minus,x,mu').^2 ./ (2* ones(numel(x),1) * SqSigma' ))) * w);

cntType = 1;

h = waitbar(0,'Please wait...');

% loop over different depth dose components  
for ixType = 1:length(Param.Types)
 
   cntComp           = 1;
   Param.Type        = Param.Types{1,ixType};
   
   % loop over the number of different Guassian components  
   for i = 1:length(vNumComp)
   
      cntEngy = 1;
      
      % loop over different energies
      for ixEnergy = vEnergyIx

         Param.energy   = machine.data(ixEnergy).energy;  
         Param.ixEnergy = ixEnergy;
         Param.NumComp  = vNumComp(cntComp);
         waitbar(ixEnergy / numel(vEnergyIx))
         
         % loop over different alpha beta ratios
         for ixBio = 1:length(Param.bioIx)
            
            switch Param.Type
                case {'Z'}

                   if ixBio > 1
                      break;
                   end
                   depths      = machine.data(ixEnergy).depths';
                   groundTruth = machine.data(ixEnergy).Z';           

                case {'alphaDose'}

                    if isstruct(machine.data(ixEnergy).Z)
                        idd = SumGauss(machine.data(ixEnergy).depths,machine.data(ixEnergy).Z.mean,...
                                                                    (machine.data(ixEnergy).Z.width).^2,...
                                                                     machine.data(ixEnergy).Z.weight);
                    else
                        idd = machine.data(ixEnergy).Z;
                    end
                    depths      = machine.data(ixEnergy).depths';
                    groundTruth = idd .* machine.data(ixEnergy).alpha(:,Param.bioIx(ixBio));

                case {'SqrtBetaDose'}

                    if isstruct(machine.data(ixEnergy).Z)
                        idd = SumGauss(machine.data(ixEnergy).depths,machine.data(ixEnergy).Z.mean,...
                                                                    (machine.data(ixEnergy).Z.width).^2,...
                                                                     machine.data(ixEnergy).Z.weight);
                    else
                        idd = machine.data(ixEnergy).Z;
                    end
                    depths      = machine.data(ixEnergy).depths';
                    groundTruth = idd .* sqrt(machine.data(ixEnergy).beta(:,Param.bioIx(ixBio)));
                   
                case {'LET'}
                    
                   if ixBio > 1
                      break;
                   end
                    groundTruth = machine.data(ixEnergy).LET';
                    depths      = machine.data(ixEnergy).depths';
                    
                    ix          = find(groundTruth>0,1,'last');
                    groundTruth = groundTruth(1:ix);
                    depths      = depths(1:ix);
                    
            end

            [w_fit, mu_fit, sigma_fit, Error] = matRad_APMfitGausComp(depths,groundTruth,Param,[],visBool);

            FitParam.(Param.Type){cntComp,cntEngy,ixBio}.weight = w_fit;
            FitParam.(Param.Type){cntComp,cntEngy,ixBio}.mean   = mu_fit;
            FitParam.(Param.Type){cntComp,cntEngy,ixBio}.width  = sigma_fit;

            mDataErr.(Param.Type)(cntComp,cntEngy,ixBio,1)   = ixEnergy;
            mDataErr.(Param.Type)(cntComp,cntEngy,ixBio,2)   = Error.rel.maxDiff;
            mDataErr.(Param.Type)(cntComp,cntEngy,ixBio,3)   = Error.rel.meanDiff;

            mDataAvg.(Param.Type)(cntComp,:,cntEngy,ixBio) = Error.rel.vDiff;

    
         end
         
         cntEngy = cntEngy + 1;
                 
      end

      cntComp   = cntComp + 1;
   end
   
   cntType = cntType + 1;
   
end
close(h)


%% plot the error
col = length(Param.Types);

try
    
   figure('Color',[1 1 1]); set(gcf, 'Position', get(0, 'Screensize'));
   defaultFontSize = 18;
   ixBio           = 1;  %alpha beta reference ratio

   for ixCol = 1:col

      for i = 1:length(vNumComp)
         vX = [machine.data(squeeze(mDataErr.(Param.Types{1,ixCol})(i,:,ixBio,1))).energy];
         subplot(col,3,(ixCol*3)-2),plot(vX,squeeze(mDataErr.(Param.Types{1,ixCol})(i,:,ixBio,2)),'LineWidth',2,'Marker','v','MarkerSize',3),hold on
         cNumComp{i} = [num2str(vNumComp(i)) '  comp. '];
      end
      xlabel('energy [MeV]','Interpreter','Latex'),ylabel('max relative difference [\%]','Interpreter','Latex'),
      title([Param.Types{1,ixCol} ' component'],'Interpreter','Latex'),grid on
      legend(cNumComp,'Interpreter','Latex','location','northwest'),set(gca,'FontSize',defaultFontSize),grid minor
      set(gca,'XLim',[min(vX) max(vX)])

      for i = 1:length(vNumComp)
          subplot(col,3,(ixCol*3)-1),plot(vX,squeeze(mDataErr.(Param.Types{1,ixCol})(i,:,ixBio,3)),'LineWidth',2,'Marker','v','MarkerSize',3),hold on
      end
      xlabel('energy [MeV]','Interpreter','Latex'),ylabel('mean relative difference [\%]','Interpreter','Latex'),grid on
      set(gca,'FontSize',defaultFontSize),grid minor,set(gca,'XLim',[min(vX) max(vX)])

      ixComp = find(vNumComp==10);
      subplot(col,3,3*ixCol),hold on
      H(1) = shadedErrorBar(linspace(0,1.6,600),squeeze(mDataAvg.(Param.Types{1,ixCol})(ixComp,:,:,ixBio))',{@mean, @(x) 2*std(x)},   '-r',0);
      H(2) = shadedErrorBar(linspace(0,1.6,600),squeeze(mDataAvg.(Param.Types{1,ixCol})(ixComp,:,:,ixBio))',{@mean, @(x) 1*std(x)},   '-g',0);
      H(3) = shadedErrorBar(linspace(0,1.6,600),squeeze(mDataAvg.(Param.Types{1,ixCol})(ixComp,:,:,ixBio))',{@mean, @(x) 0.5*std(x)}, '-b',0);
      grid on,grid minor
      legend([H(3).mainLine, H.patch], '\mu','2\sigma','\sigma','0.5\sigma');
      xlabel('relative depth','Interpreter','Latex'),title(['approx error using ' num2str(vNumComp(ixComp)) ' components' ],'Interpreter','Latex');
      ylabel('relative difference  [\%]','Interpreter','Latex'),set(gca,'FontSize',defaultFontSize);

   end
   
catch
    
end

save('apmFIT.mat')

%% save data to new machine file
load('apmFIT.mat')

SumGauss = @(X,MU,SqSigma,W) ((1./sqrt(2*pi*ones(numel(X),1) * SqSigma') .* ...
                              exp(-bsxfun(@minus,X,MU').^2 ./ (2* ones(numel(X),1) * SqSigma' ))) * W);
SG       =  @(qX,qW,qMu,qSigma)((qW/(sqrt(2*pi*qSigma^2))).*exp(-((qX-qMu).^2)./(2*qSigma^2)));  

% copy old file
machineRef   = machine;
ixComp       = [];
col          = length(Param.Types);
if isempty(ixComp)
    ixComp = 1;
end

for ixCol = 1:col
   
   if isfield(machine.data,Param.Types{1,ixCol})
      machine.data = rmfield(machine.data,Param.Types{1,ixCol});
   end
   
   if isequal(Param.Types{1,ixCol},'alphaDose') 
       machine.data = rmfield(machine.data,'alphaX');
       machine.data = rmfield(machine.data,'alpha');
       machine.data = rmfield(machine.data,'alphaBetaRatio');
   end
   
   if isequal(Param.Types{1,ixCol},'SqrtBetaDose')
       machine.data = rmfield(machine.data,'betaX');
       machine.data = rmfield(machine.data,'beta');
   end
   
   CntEnergy = 1;
   for i = vEnergyIx
      
      CntBio = 1;
      
      for ixBio = Param.bioIx
  
          if isequal(Param.Types{1,ixCol},'Z')
             idd  = machineRef.data(i).Z;
             if CntBio > 1
                break;
             end
          elseif isequal(Param.Types{1,ixCol},'alphaDose')
             idd  = machineRef.data(i).Z .* machineRef.data(i).alpha(:,ixBio);
          elseif isequal(Param.Types{1,ixCol},'SqrtBetaDose')
             idd  = machineRef.data(i).Z .* sqrt(machineRef.data(i).beta(:,ixBio));
          elseif  isequal(Param.Types{1,ixCol},'LET')
             idd  = machineRef.data(i).LET; 
             if ixBio > 1
                break;
             end
          end

          machine.data(i).(Param.Types{1,ixCol})(CntBio).weight      = FitParam.(Param.Types{1,ixCol}){ixComp,CntEnergy,CntBio}.weight';
          machine.data(i).(Param.Types{1,ixCol})(CntBio).mean        = FitParam.(Param.Types{1,ixCol}){ixComp,CntEnergy,CntBio}.mean';
          machine.data(i).(Param.Types{1,ixCol})(CntBio).width       = FitParam.(Param.Types{1,ixCol}){ixComp,CntEnergy,CntBio}.width';
          
          
          machine.data(i).(Param.Types{1,ixCol})(CntBio).profileORG  = idd;
          machine.data(i).(Param.Types{1,ixCol})(CntBio).profileAPM  = SumGauss(machine.data(i).depths,...
                                                             FitParam.(Param.Types{1,ixCol}){ixComp,CntEnergy,CntBio}.mean',...
                                                             FitParam.(Param.Types{1,ixCol}){ixComp,CntEnergy,CntBio}.width'.^2,...
                                                             FitParam.(Param.Types{1,ixCol}){ixComp,CntEnergy,CntBio}.weight');
          
          
          if isequal(Param.Types{1,ixCol},'alphaDose')
             machine.data(i).alphaX(1,CntBio) = machineRef.data(i).alphaX(1,Param.bioIx(CntBio));
             machine.data(i).alpha(:,CntBio)  = machineRef.data(i).alpha(:,Param.bioIx(CntBio));
             machine.data(i).alphaBetaRatio(1,CntBio) = machineRef.data(i).alphaBetaRatio(1,Param.bioIx(CntBio));
          end
          
          if isequal(Param.Types{1,ixCol},'SqrtBetaDose')
             machine.data(i).betaX(1,CntBio) = machineRef.data(i).betaX(1,Param.bioIx(CntBio));
             machine.data(i).beta(:,CntBio)  = machineRef.data(i).beta(:,Param.bioIx(CntBio));
          end
          
          CntBio = CntBio + 1;  
      end
      
      CntEnergy = CntEnergy + 1;      
      
   end
   
end


%% save result
save([machine.meta.radiationMode '_' machine.meta.machine 'APM'],'machine')

%%
if visAllBool
   
   % plot single energies
   figure('Color',[1 1 1]); cSpecs = colorspecs;
   ix   = 11;  Comp = 'Z'; ixBio = 1;
   xx = 0:0.2:machine.data(ix).depths(end);
   yy = zeros(length(xx),1);

   for i = 1:numel(machine.data(ix).(Comp).weight)   
       u = SG(xx,machine.data(ix).(Comp)(ixBio).weight(i),machine.data(ix).(Comp)(ixBio).mean(i),machine.data(ix).(Comp)(ixBio).width(i));
       plot(xx,u,'LineWidth',2,'color',cSpecs.dre),hold on
       yy = yy + u';
   end
   plot(xx,yy,'LineWidth',3,'color',cSpecs.dkfzdB)
   plot(machine.data(ix).depths,machine.data(ix).(Comp)(1).profileORG,'--','LineWidth',3,'color',cSpecs.dkfzlB), 

   xlabel('x [mm]','Interpreter','Latex'),ylabel('rel. dose','Interpreter','Latex')

   set(gca,'xlim',[0 xx(end)])
   title(['func. approx. ' Comp ', E = ' num2str(machine.data(ix).energy) ' MeV'],'Interpreter','Latex')
   set(gca,'FontSize',18);

   % plot fits of all energies
   figure('Color',[1 1 1]);
   for i=1:1:length(machine.data)
       cla
       x = machine.data(i).depths;
       y1 = zeros(length(x),1);
       y = SumGauss(x,machine.data(i).Z.mean,machine.data(i).Z.width.^2,machine.data(i).Z.weight);

       for j = 1:length(vNumComp(ixComp))
          y1 = y1 + machine.data(i).Z.weight(j) * normpdf(x,machine.data(i).Z.mean(j),machine.data(i).Z.width(j));
       end

       plot(x,y,'r','LineWidth',3),hold on
       plot(x,machineRef.data(i).Z,'b','LineWidth',2);
       legend({'fitted','ground truth'}),title('press a button')
       waitforbuttonpress 
   end


   figure('Color',[1 1 1]);
   for i=1:1:length(machine.data)
       cla
       x = (0:0.2:600)';
       y1 = zeros(length(x),1);

       y = SumGauss(x,machine.data(i).Z.mean,machine.data(i).Z.width.^2,machine.data(i).Z.weight);
       for j = 1:(vNumComp(ixComp))
          ysg = machine.data(i).Z.weight(j) * normpdf(x,machine.data(i).Z.mean(j),machine.data(i).Z.width(j));
          plot(x,ysg),hold on
          y1  = y1 + ysg;
       end
       plot(x,y,'k','LineWidth',2),hold on,title('press a button')
       waitforbuttonpress 
   end


end
