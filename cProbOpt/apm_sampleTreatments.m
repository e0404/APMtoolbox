function doseSampleS = apm_sampleTreatments(nS_S, nFrac, mu, sigma, ucm, xSamp, w, showWaitbar, showScatterPlots)
% Generates nS_S treatments of nFrac fractions each, with systematic and
% random uncertainties in the positions of the spots of weight w, 
% evaluated at the voxels in positions xSamp

if nargin < 8
    showWaitbar = true;
end

if nargin < 9
    showScatterPlots = false;
end

% Sample the systematics in the positions of the spots (same for all fractions)
S_S = mb_mgd(nS_S,mu,ucm.covSys)'; 

% Create the matrix that will store the sampled treatments
doseSampleS = zeros(numel(xSamp),nS_S);

for i = 1:nS_S
    % Sample the random fluctuations in the positions of the spots (that are
    % already shifted from the systematics)
    S_R = mb_mgd(nFrac,S_S(:,i)',ucm.covRand)'; 
    
    doseSampleF = zeros(numel(xSamp), nFrac);

    for f = 1:nFrac
        shifted_Spots = struct('mu',num2cell(S_R(:,f)),'sigma', num2cell(sigma)); 
        
        % Calculate dose at voxels in xSamp
        doseSampleF_ij = apm_calcDoseInfluenceLateral(xSamp,shifted_Spots);
        doseSampleF(:,f) = apm_calcDose(doseSampleF_ij,w);
    end
    doseSampleS(:,i) = sum(doseSampleF,2)/nFrac;
    
    % Plot the accumulated treatment dose sample
    if (showScatterPlots)
        plot(xSamp,doseSampleS(:,i),'.','Color',[0.6 0.6 1],'MarkerSize',2);
        drawnow;
    end
    
    %plot(xSamp,std(yR),'.','Color',[1 0.6 1],'MarkerSize',2)
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
