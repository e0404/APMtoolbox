function [expDvh,stdDvh,covDvh] = apm_DVHprob(expDose,covDose,nBins,dMax,method)

if nargin < 3 || isempty(nBins)
    nBins = 100;
end

if nargin < 4 || isempty(dMax)
    dMax = 1.25*max(expDose);
end

if nargin < 5
    method = 'int_gauss';
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
    case 'int_gauss'
        [expDvh(2,:),stdDvh(2,:)] = probDVHGaussExpStd(expDose,covDose,dAxis);
        if nargout > 2
            covDvh = calcCovDvhGauss(expDose,covDose,dAxis,expDvh(2,:));
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
            %P2D = arrayfun(@(dParam) mexBVNcdf(-[dParam dParam],-[expDose(i) expDose(l)],[covDose(i,i) covDose(i,l); covDose(l,i) covDose(l,l)]),dAxis);
            
            if any(isnan(P2D))
                warn('CDF (2D) evaluates to NaN');
                P2D(isnan(P2D)) = 0;
            end
            
            biCum = biCum + P2D;
        end
        matRad_progress(i,numel(expDose));
    end
    
    expDvh = diagCum ./ numel(expDose);
    stdDvh = 1/numel(expDose)^2 * (diagCum + 2*biCum);

    stdDvh = sqrt(stdDvh - expDvh.^2);
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
