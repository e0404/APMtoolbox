function [doseSampleS, nS_S] = apm_runSampleValidation(w, axDVH, expDose, stdDose, nSamplesTotal, nFrac, mu, sigma, ucm, xLowRes, ixLowRes, xStar, axProfile, t, probDvhMethod, nDvhBins, vois, showWaitbar, showScatterPlots, showAllSampDVHs, showAvgSampDVHs)
    % Simulate a number of treatments with fractionation.
    % In the profile plot, plot the mean dose +- sdev of the mean,
    % and stdev +- stdev of the stdev, and for a specific xStar, 
    % plot the distribution of the dose.
    % Separately, plot the nominal DVHs of each sample.
    
    if nargin<18
        showWaitbar=true;
    end
    if nargin<19
        showScatterPlots=false;
    end
    % Compute number of treatments
    nS_S = floor(nSamplesTotal/nFrac);

    % Sample the treatments
    doseSampleS = apm_sampleTreatments(nS_S, nFrac, mu, sigma, ucm, xLowRes, w, showWaitbar, showScatterPlots);
    doseSampleS_xStar = apm_sampleTreatments(nS_S, nFrac, mu, sigma, ucm, xStar, w, showWaitbar, showScatterPlots);

    % Compute statistics
    doseSampleS_mean = mean(doseSampleS,2);
    doseSampleS_std  = std(doseSampleS,[],2);
    %doseSampleS_cov = cov(doseSampleS');
    doseSampleS_mean_stdErr = doseSampleS_std / sqrt(nS_S);
    doseSampleS_std_stdErr  = doseSampleS_std * apm_calcSampledStdAccRel(nS_S);

    doseSampleS_xStar_mean = mean(doseSampleS_xStar,2);
    doseSampleS_xStar_std  = std(doseSampleS_xStar,[],2);
    
    % Plot mean and stdev dose as a function of x 
    errorbar(axProfile,xLowRes,doseSampleS_mean,doseSampleS_mean_stdErr,'rs','MarkerSize',2);
    errorbar(axProfile,xLowRes,doseSampleS_std,doseSampleS_std_stdErr,'ms','MarkerSize',2);
    %plot(xLowRes,kurtosis(yS),'y*','MarkerSize',10);
    %plot(xLowRes,skewness(yS),'k*','MarkerSize',10);
    

    % Plot normalized histogram of dose at xStar
    axProfile2 = axes(t);
    [counts, centers] = hist(doseSampleS_xStar,floor(sqrt(nS_S)));
    binWidths = diff(centers);
    counts = (counts/binWidths(1))/nS_S;
    % Adjust bar heights to be contained within plot limits
    %if ((max(counts)+xStar)>max(vois(1).xU, vois(2).xU))
    %    counts = counts*(max(vois(1).xU, vois(2).xU))*0.75/(max(counts)+xStar); 
    %end
    %counts = counts*((max(vois(1).xU, vois(2).xU) - min(vois(1).xL, vois(2).xL)))*0.15/(max(counts));
    barh(axProfile2,centers, counts,'EdgeColor','none','BarWidth',1,'FaceColor',.8*[1 1 1]);%counts+xStar, 'BaseValue', xStar);
    plot(axProfile,[xStar xStar],[0 1.2],'--','Color','k','MarkerFaceColor','k');
    
    % Plot modelled distribution of dose at xStar  
    axProfile3 = axes(t);
    maxPdf = apm_plotModelledDosePdf(axProfile3, probDvhMethod, doseSampleS_xStar_mean, doseSampleS_xStar_std);
    
    % Make the histogram and modelled distribution occupy the left 15% of the plot
    xlim(axProfile3, [0, max(real(max(counts)/0.15), real(maxPdf/0.15))]); 
    axProfile2.XAxisLocation = 'top';
    axProfile3.XAxisLocation = 'top';
    xlabel(axProfile3, 'rel. counts');
    axProfile2.Color = 'none';
    axProfile3.Color = 'none';
    axProfile.Box = 'off';
    axProfile2.Box = 'off';
    axProfile3.Box = 'off';
    linkaxes([axProfile axProfile2 axProfile3], 'y');
    linkaxes([axProfile2 axProfile3], 'x');

    % Compute nominal DVHs for the samples
    disp(['Computing nominal DVHs for samples']); 
    sampDVH_x = zeros(1, nDvhBins);
    sampDVH_y_v1 = zeros(1, nDvhBins);
    sampDVH_y_v2 = zeros(1, nDvhBins);
    for i = 1:nS_S    
        for v = 1:numel(vois)
            % Get the indices of the voxels in ixLowRes that correspond to
            % those of vois(v)
            isxLowResInVoi = ismember(ixLowRes,find(vois(v).ix==1));
            idxs = find(isxLowResInVoi==1);
            % Compute the nominal DVHs
            sampDVH_iv = apm_DVH(doseSampleS(idxs,i),nDvhBins,1.1);
            % Get the values of the dose (same for all vois and samples)
            sampDVH_x = sampDVH_iv(1,:); 
            if v==1
                sampDVH_y_v1 = vertcat(sampDVH_y_v1, sampDVH_iv(2,:));
            end
            if v==2
                sampDVH_y_v2 = vertcat(sampDVH_y_v2, sampDVH_iv(2,:));
            end
            %plot(axSamples,sampDVH_iv(1,:),sampDVH_iv(2,:),'LineWidth',1);
        end
    end
    sampDVH_y_v1 = sampDVH_y_v1(2:(nS_S+1),:);
    sampDVH_y_v2 = sampDVH_y_v2(2:(nS_S+1),:);
    
    % Compute statistics on the samples
    sampDVH_v1_mean = mean(sampDVH_y_v1,1);
    sampDVH_v1_std  = std(sampDVH_y_v1,[],1);
    sampDVH_v2_mean = mean(sampDVH_y_v2,1);
    sampDVH_v2_std  = std(sampDVH_y_v2,[],1);

    % Plot the individual sampled DVHs
    if showAllSampDVHs
        figure;
        axSamples = axes();
        hold(axSamples, 'on');
        for i = 1:nS_S    
                plot(axSamples,sampDVH_y_v1(i,:),sampDVH_y_v2(i,:),'LineWidth',1);
        end
        box(axSamples,'on');
        grid(axSamples,'on');
        ylim(axSamples,[0 1]);
        xlim(axSamples,[0 1.1]);
        xlabel(axSamples,'rel. dose');
        ylabel(axSamples,'rel. volume');
        title('Sample DHVs');
    end

    % Plot the average sampled DVHs
    if showAvgSampDVHs
        plot(axDVH, sampDVH_x, sampDVH_v1_mean, 'ks', 'MarkerSize',0.5);
        plot(axDVH, sampDVH_x, sampDVH_v2_mean, 'ks', 'MarkerSize',0.5);
        fill(axDVH, [sampDVH_x'; flipud(sampDVH_x')], [(sampDVH_v1_mean - sampDVH_v1_std)'; flipud((sampDVH_v1_mean + sampDVH_v1_std)')],'k','FaceAlpha',0.1,'LineStyle','none');
        fill(axDVH, [sampDVH_x'; flipud(sampDVH_x')], [(sampDVH_v2_mean - sampDVH_v2_std)'; flipud((sampDVH_v2_mean + sampDVH_v2_std)')],'k','FaceAlpha',0.1,'LineStyle','none');
    end
    
    %Evaluate accuracy of sampled vs analytical dose(x)
    dose_mean_atSamples = expDose(ixLowRes);
    dose_std_atSamples = stdDose(ixLowRes);
    
    rmse_mean = sqrt(mean((dose_mean_atSamples - doseSampleS_mean).^2));
    rmse_std  = sqrt(mean((dose_std_atSamples - doseSampleS_std).^2));
    
    disp(['Mean RMSE: ' num2str(rmse_mean)]);
    disp(['Std RMSE: ' num2str(rmse_std)]);
end