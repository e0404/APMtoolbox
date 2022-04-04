function [expDvh,stdDvh,covDvh] = apm_DVHprob(expDose,covDose,nBins,dMax,method,copula_model,copula_marginals,dmin_shiftedbeta,dmax_shiftedbeta)

if nargin < 3 || isempty(nBins)
    nBins = 100;
end

if nargin < 4 || isempty(dMax)
    dMax = 1.25*max(expDose);
end

if nargin < 5
    method = 'int_gauss';
end

% Define parameters of the Copula
if nargin < 6
    copula_model = 'gauss';
end

if nargin < 7
    copula_marginals = cell(1,numel(expDose));
    copula_marginals(:) = {'gauss'};
end

% Define support of the shifted beta
if nargin < 8
    dmin_shiftedbeta= zeros(1,numel(expDose));
end

if nargin < 9
    dmax_shiftedbeta = ones(1, numel(expDose));
end

dAxis = linspace(0,dMax,nBins);
expDvh(1,:) = dAxis;
stdDvh(1,:) = dAxis;

fullVoxNum = numel(expDose);

%{
% remove the zero variance voxels
varDose = diag(covDose);
ix = find(varDose);
covDose = covDose(ix,ix);
expDose = expDose(ix);
%}

% make variance non-zero
varDose = diag(covDose);
addDiag = zeros(size(varDose));
addDiag(varDose == 0) = 1e-9;
covDose = covDose + diag(addDiag);

reducedVoxNum = numel(expDose);

switch method
    case 'copula'
        [expDvh(2,:),stdDvh(2,:)] = probDVHCopulaExpStd(expDose,covDose,dAxis,copula_model,copula_marginals,dmin_shiftedbeta,dmax_shiftedbeta);
        if nargout > 2
            covDvh = calcCovDvhCopula(expDose,covDose,dAxis,expDvh(2,:),copula_model,copula_marginals,dmin_shiftedbeta,dmax_shiftedbeta);
        end
    case 'int_gauss'
        [expDvh(2,:),stdDvh(2,:)] = probDVHGaussExpStd(expDose,covDose,dAxis);
        if nargout > 2
            covDvh = calcCovDvhGauss(expDose,covDose,dAxis,expDvh(2,:));
        end

    case 'int_loggauss'
        [expDvh(2,:),stdDvh(2,:)] = probDVHLnGaussExpStd(expDose,covDose,dAxis);
        if nargout > 2
            covDvh = calcCovDvhLnGauss(expDose,covDose,dAxis,expDvh(2,:));
        end
    
    case 'int_gauss_cuda'
        disp('Setting up kernel');
        tic
        

        %single precision
        %{
        cuKernel = parallel.gpu.CUDAKernel('C:\dev\APM\CUDA\QI\DVH\x64\Release\dvh.ptx','float* , float*, const float*, const float*, const float*, int, int','dvhRawMoments');
        cuKernel.ThreadBlockSize = [16 16 3];
        firstRawMomentGPU = gpuArray(zeros(nBins,1,'single'));
        secondRawMomentGPU = gpuArray(zeros(nBins,1,'single'));
        dAxisGPU = gpuArray(single(transpose(dAxis(:))));
        expDoseGPU = gpuArray(single(expDose));
        covDoseGPU = gpuArray(single(reshape(covDose,[nVox*nVox 1])));
        %}
        
        %double precision
        cuKernelExp = parallel.gpu.CUDAKernel('C:\dev\APM\CUDA\QI\DVH\x64\Release\dvh.ptx','double* , const double*, const double*, const double*, int, int','dvhExp');
        cuKernelExp.ThreadBlockSize = [128 8 1];
        
        %cuKernelVar = parallel.gpu.CUDAKernel('C:\dev\APM\CUDA\QI\DVH\x64\Release\dvh.ptx','double* , const double*, const double*, const double*, const double*, int, int','dvhVar');
        %cuKernelVar.ThreadBlockSize = [12 12 4]; 
        
        cuKernelVarSingle = parallel.gpu.CUDAKernel('C:\dev\APM\CUDA\QI\DVH\x64\Release\dvh.ptx','double* , const double*, const double*, const double*, const double*, int, int','dvhVarSingleBin');
        cuKernelVarSingle.ThreadBlockSize = [24 24 1]; 
        
        %cuKernelCov = parallel.gpu.CUDAKernel('C:\dev\APM\CUDA\QI\DVH\x64\Release\dvh.ptx','double* , const double*, const double*, const double*, const double*, int, int','dvhCov');
        %cuKernelCov.ThreadBlockSize = [12 12 4]; 
        cuKernelCovSingle = parallel.gpu.CUDAKernel('C:\dev\APM\CUDA\QI\DVH\x64\Release\dvh.ptx','double* , const double*, const double*, const double*, const double*, int, int, int, int','dvhCovTwoBins');
        cuKernelCovSingle.ThreadBlockSize = [24 24 1]; 
        
        expDvhGPU = gpuArray(zeros(nBins,1,'double'));
        varDvhGPU = gpuArray(zeros(nBins,1,'double'));
        %covDvhGPU = gpuArray(zeros(nBins,nBins,'double')); 
        dAxisGPU = gpuArray(double(transpose(dAxis(:))));
        expDoseGPU = gpuArray(double(expDose'));
        %covDoseGPU = gpuArray(double(reshape(covDose,[nVox*nVox 1])));
        covDoseGPU = gpuArray(double(covDose));
        varDoseGPU = gpuArray(double(diag(covDose)));
        
        
        cuKernelExp.GridSize = ceil([reducedVoxNum nBins 1] ./ cuKernelExp.ThreadBlockSize);        
        %cuKernelVar.GridSize = ceil([nVox nVox nBins] ./ cuKernelVar.ThreadBlockSize);    
        cuKernelVarSingle.GridSize = ceil([reducedVoxNum reducedVoxNum 1] ./ cuKernelVarSingle.ThreadBlockSize);
        %cuKernelCov.GridSize = ceil([nVox nVox nBins^2] ./ cuKernelCov.ThreadBlockSize);        
        cuKernelCovSingle.GridSize = ceil([reducedVoxNum reducedVoxNum 1] ./ cuKernelCovSingle.ThreadBlockSize);        
        
        
        disp(['    CUDA Grid size for E[DVH]:   ' num2str(cuKernelExp.GridSize)]);
        %disp(['    CUDA Grid size for Var[DVH]: ' num2str(cuKernelVar.GridSize)]);
        disp(['    CUDA Grid size for Var[DVH]: ' num2str(cuKernelVarSingle.GridSize)]);
        %disp(['    CUDA Grid size for Cov[DVH]: ' num2str(cuKernelCov.GridSize)]);
        disp(['    CUDA Grid size for Cov[DVH]: ' num2str(cuKernelCovSingle.GridSize)]);
        
        disp('Running kernels...');
        %[firstRawMomentGPU,secondRawMomentGPU] = feval(cuKernelExp,firstRawMomentGPU,secondRawMomentGPU,dAxisGPU,expDoseGPU,covDoseGPU,nVox,nBins);
        tic
        expDvhGPU = feval(cuKernelExp,expDvhGPU,expDoseGPU,varDoseGPU,dAxisGPU,reducedVoxNum,nBins);
        wait(gpuDevice)
        toc
        
        if nargout > 2
            covDvhGPU = gpuArray(zeros(nBins,nBins,'double')); 
            tic
            %covDvhGPU = feval(cuKernelCov,covDvhGPU,expDoseGPU,covDoseGPU,dAxisGPU,expDvhGPU,nVox,nBins);
            c=1;
            for p = 1:nBins
                for q=1:nBins
                    covDvhGPU = feval(cuKernelCovSingle,covDvhGPU,expDoseGPU,covDoseGPU,dAxisGPU,expDvhGPU,reducedVoxNum,nBins,p-1,q-1);
                    wait(gpuDevice);
                    matRad_progress(c,nBins^2);
                    c=c+1;
                end
            end
            wait(gpuDevice)
            toc
            varDvhGPU = diag(covDvhGPU);
            covDvh = gather(covDvhGPU);          
        else
            tic
            %varDvhGPU = feval(cuKernelVar,varDvhGPU,expDoseGPU,covDoseGPU,dAxisGPU,expDvhGPU,nVox,nBins);            
            for dIx = 1:nBins
                varDvhGPU = feval(cuKernelVarSingle,varDvhGPU,expDoseGPU,covDoseGPU,dAxisGPU,expDvhGPU,reducedVoxNum,dIx-1);
                wait(gpuDevice);
                matRad_progress(dIx,nBins);
            end
            wait(gpuDevice);
            toc
        end
                    
        expDvh(2,:) = transpose(gather(expDvhGPU));
        stdDvh(2,:) = sqrt(transpose(gather(varDvhGPU)));
        
    otherwise
        error(['Method ''' method ''' for computing probabilistic DVH not known']);
        expDvh(2,:) = NaN*zeros(1,nBins);
        stdDvh(2,:) = NaN*zeros(1,nBins);
        covDvh = NaN*zeros(nBins);        
end

expDvh(2,:) = expDvh(2,:)*reducedVoxNum/fullVoxNum;

end

function [expDvh,stdDvh] = probDVHCopulaExpStd(expDose,covDose,dAxis,copula_model,copula_marginals,dmin_shiftedbeta,dmax_shiftedbeta)    
    nBins = numel(dAxis);
    biCum = zeros(1,nBins);
    diagCum = zeros(1,nBins);

    for i = 1:numel(expDose)
        mean_i = expDose(i);
        sig_i = sqrt(covDose(i,i));        
        if ~isreal(sig_i) || isnan(sig_i)
            continue;
        end
        P1 = apm_calcCdfDose(dAxis, mean_i*ones(numel(1,nBins)), sig_i*ones(numel(1,nBins)), copula_marginals(i), 1,dmin_shiftedbeta(i),dmax_shiftedbeta(i));
        
        
        if any(isnan(P1))
            warning('CDF evaluates to NaN');
            P1(isnan(P1)) = 0;
        end
        diagCum = diagCum + P1;
        
        for l = i+1:numel(expDose)
            mean_l = expDose(l);
            sig_l = sqrt(covDose(l,l));
            if ~isreal(sig_l) || isnan(sig_l)
                continue;
            end
            r = (covDose(i,l) / (sig_i*sig_l)) * ones(1,nBins);
            u_i = apm_calcCdfDose(dAxis, mean_i, sig_i, copula_marginals(i),0,dmin_shiftedbeta(i),dmax_shiftedbeta(i));
            u_l = apm_calcCdfDose(dAxis, mean_l, sig_l, copula_marginals(l),0,dmin_shiftedbeta(l),dmax_shiftedbeta(l));
            P2D = apm_calcCopula(u_i, u_l, r, copula_model) + 1 - u_i - u_l;
            
            if any(isnan(P2D))
                warning('CDF (2D) evaluates to NaN');
                P2D(isnan(P2D)) = 0;
            end
            
            biCum = biCum + P2D;
        end
        matRad_progress(i,numel(expDose));
    end

    % biCum has only the terms from l=i+1 to numel(expDose), 
    % then we need to add the terms l=i (diagCum) and l from 0 to i-1 (biCum).
    biCum = diagCum + 2*biCum;
    
    expDvh = diagCum ./ numel(expDose);

    stdDvh = sqrt( 1/numel(expDose)^2 * (biCum) - expDvh.^2);
end

function covDvh = calcCovDvhCopula(expDose,covDose,dAxis,expDvh,copula_model,copula_marginals,dmin_shiftedbeta,dmax_shiftedbeta)
    nBins = numel(dAxis);
    covDvh = zeros(nBins);
    
    for p = 1:nBins
        for q = 1:nBins
            s = 0;
            for i = 1:numel(expDose)
                mean_i = expDose(i);
                sig_i = sqrt(covDose(i,i));
                for l = 1:numel(expDose)
                    mean_l = expDose(l);
                    sig_l = sqrt(covDose(l,l));
                    r = covDose(i,l) / (sig_i*sig_l);
                    u_i = apm_calcCdfDose(dAxis(p), mean_i, sig_i, copula_marginals(i),0,dmin_shiftedbeta(i),dmax_shiftedbeta(i));
                    u_l = apm_calcCdfDose(dAxis(q), mean_l, sig_l, copula_marginals(l),0,dmin_shiftedbeta(l),dmax_shiftedbeta(l));
                    s = s + apm_calcCopula(u_i, u_l, r, copula_model) + 1 - u_i - u_l;
                end
            end
            
            covDvh(p,q) = s;
        end
        matRad_progress(p,nBins);
    end
    covDvh = 1/numel(expDose)^2 * covDvh;
    covDvh = covDvh - transpose(expDvh)*expDvh;
end


function [expDvh,stdDvh] = probDVHGaussExpStd(expDose,covDose,dAxis)
    nBins = numel(dAxis);
    biCum = zeros(1,nBins);
    diagCum = zeros(1,nBins);
    
    for i = 1:numel(expDose)
        sig_i = sqrt(covDose(i,i));
        
        if ~isreal(sig_i) || isnan(sig_i)
            continue;
        end
        
        P1 = 1 - apm_normcdf(dAxis,expDose(i)*ones(numel(1,nBins)),sig_i*ones(numel(1,nBins)));
        
        
        if any(isnan(P1))
            warning('CDF evaluates to NaN');
            P1(isnan(P1)) = 0;
        end
        diagCum = diagCum + P1;
        
        for l = i+1:numel(expDose)
            sig_l = sqrt(covDose(l,l));
            
            if ~isreal(sig_l) || isnan(sig_l)
                continue;
            end
            
            r = covDose(i,l) / (sig_i*sig_l);
            P2D = arrayfun(@(dParam) bvn((dParam - expDose(i)) / sig_i,Inf,(dParam - expDose(l)) / sig_l,Inf,r),dAxis);
            
            if any(isnan(P2D))
                warn('CDF (2D) evaluates to NaN');
                P2D(isnan(P2D)) = 0;
            end
            
            biCum = biCum + P2D;
        end
        matRad_progress(i,numel(expDose));
    end

    % biCum has only the terms from l=i+1 to numel(expDose), 
    % then we need to add the terms l=i (diagCum) and l from 0 to i-1 (biCum).
    biCum = diagCum + 2*biCum;
    
    expDvh = diagCum ./ numel(expDose);

    stdDvh = sqrt( 1/numel(expDose)^2 * (biCum) - expDvh.^2);
end

function covDvh = calcCovDvhGauss(expDose,covDose,dAxis,expDvh)
    nBins = numel(dAxis);
    covDvh = zeros(nBins);
    
    for p = 1:nBins
        for q = 1:nBins
            s = 0;
            for i = 1:numel(expDose)
                sig_i = sqrt(covDose(i,i));
                for l = 1:numel(expDose)
                    sig_l = sqrt(covDose(l,l));
                    r = covDose(i,l) / (sig_i*sig_l);
                    s = s + bvn((dAxis(p) - expDose(i)) / sig_i,Inf,(dAxis(q) - expDose(l)) / sig_l,Inf,r);
                end
            end
            
            covDvh(p,q) = s;
        end
        matRad_progress(p,nBins);
    end
    covDvh = 1/numel(expDose)^2 * covDvh;
    covDvh = covDvh - transpose(expDvh)*expDvh;
end


function [expDvh,stdDvh] = probDVHLnGaussExpStd(expDose,covDose,dAxis)
    nBins = numel(dAxis);
    biCum = zeros(1,nBins);
    diagCum = zeros(1,nBins);
    
    %Compute parameters of the lognormal
    [lnExpDose, lnCovDose] = apm_transformMeanCovarianceToLogNormalParameters(expDose,covDose);

    for i = 1:numel(lnExpDose) %for each voxel
        lnSig_i = sqrt(lnCovDose(i,i));
        
        if ~isreal(lnSig_i) || isnan(lnSig_i)
            continue;
        end
        
        P1 = 1 - apm_normcdf(log(dAxis),lnExpDose(i)*ones(numel(1,nBins)),lnSig_i*ones(numel(1,nBins)));
        
        if any(isnan(P1))
            warning('CDF evaluates to NaN');
            P1(isnan(P1)) = 0;
        end
        diagCum = diagCum + P1;
        
        for l = i+1:numel(lnExpDose)
            lnSig_l = sqrt(lnCovDose(l,l));
            
            if ~isreal(lnSig_l) || isnan(lnSig_l)
                continue;
            end
            
            r = lnCovDose(i,l) / (lnSig_i*lnSig_l);
            P2D = arrayfun(@(dParam) bvn((log(dParam) - lnExpDose(i)) / lnSig_i,Inf,(log(dParam) - lnExpDose(l)) / lnSig_l,Inf,r), dAxis);
            
            if any(isnan(P2D))
                warn('CDF (2D) evaluates to NaN');
                P2D(isnan(P2D)) = 0;
            end
            
            biCum = biCum + P2D;
        end
        matRad_progress(i,numel(expDose));
    end

    % biCum has only the terms from l=i+1 to numel(expDose), 
    % then we need to add the terms l=i (diagCum) and l from 0 to i-1 (biCum).
    biCum = diagCum + 2*biCum;
    
    expDvh = diagCum ./ numel(lnExpDose);

    stdDvh = sqrt( 1/numel(lnExpDose)^2 * (biCum) - expDvh.^2);
end

function covDvh = calcCovDvhLnGauss(expDose,covDose,dAxis,expDvh)
    nBins = numel(dAxis);
    covDvh = zeros(nBins);

    %Compute parameters of the lognormal
    [lnExpDose, lnCovDose] = apm_transformMeanCovarianceToLogNormalParameters(expDose,covDose); 
    
    for p = 1:nBins
        for q = 1:nBins
            s = 0;
            for i = 1:numel(lnExpDose)
                lnSig_i = sqrt(lnCovDose(i,i));
                for l = 1:numel(lnExpDose)
                    lnSig_l = sqrt(lnCovDose(l,l));
                    r = lnCovDose(i,l) / (lnSig_i*lnSig_l);
                    s = s + bvn((log(dAxis(p)) - lnExpDose(i)) / lnSig_i,Inf,(log(dAxis(q)) - lnExpDose(l)) / lnSig_l,inferiorto,r);
                end
            end
            
            covDvh(p,q) = s;
        end
        matRad_progress(p,nBins);
    end
    covDvh = 1/numel(lnExpDose)^2 * covDvh;
    covDvh = covDvh - transpose(expDvh)*expDvh;
end
