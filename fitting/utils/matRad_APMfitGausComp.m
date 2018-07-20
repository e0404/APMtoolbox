function [w_fit, mu_fit, sigma_fit, Error] = matRad_APMfitGausComp(qX,qY,Param,wInit,visBool)

% fitting function that fits Gaussian components to the profile defined by
% qX and qY. This function allows to loop over multiple optimizers.
% Further, matRad_APMfitGetIniParameters needs probably some manual
% adjustment of the initial fit parameter

SG       = @(qX,qW,qMu,qSigma)((qW/(sqrt(2*pi*qSigma^2))).*exp(-((qX-qMu).^2)./(2*qSigma^2)));
SumGauss = @(x,mu,SqSigma,w) ((1./sqrt(2*pi*ones(numel(x),1) * SqSigma') .* ...
    exp(-bsxfun(@minus,x,mu').^2 ./ (2* ones(numel(x),1) * SqSigma' ))) * w);

% get a finer resolution
qX_long = qX(1):0.25:qX(end);
qY_long = interp1(qX,qY,qX_long);

if isempty(wInit)
    [w_ini, mu_ini, sigma_ini] = matRad_APMfitGetIniParameters(qX,qY,Param);
    wInit    = [w_ini' mu_ini' sigma_ini'];
else
    [w_ini, mu_ini, sigma_ini] = splitOptResult(wInit);
    w_ini     = w_ini';
    mu_ini    = mu_ini';
    sigma_ini = sigma_ini';
end

NumParam = length(wInit);
NumOpti  = length(Param.Optimizer);
colors   = colorspecs;

vLowerBound                                    = ones(1,NumParam);       % Lower bound on on the weights and sigma.
vLowerBound(1+Param.NumComp:2*Param.NumComp) = -200;                     % Lower bound on the mean.
vUpperBound                                    = 1e4 * ones(1,NumParam); % Upper bound on  weight sigma and mean.

if visBool
    figure('Color',[1 1 1]);
end

for ixOpt = 1:NumOpti
    
    currOpt = Param.Optimizer{1,ixOpt};
    
    if visBool
        subplot(NumOpti,2,(ixOpt*2)-1),plot(qX_long,qY_long,'color',colors.dre,'LineWidth',2),grid on, grid minor, hold on
        for i = 1:numel(w_ini)
            subplot(NumOpti,2,(ixOpt*2)-1),plot(qX_long,SG(qX_long,w_ini(i),mu_ini(i),sigma_ini(i)),'color',colors.dre),hold on
        end
        set(gca,'Xlim',[min(qX_long) max(qX_long)])
        
        subplot(NumOpti,2,(ixOpt*2)-1),plot(qX_long,SumGauss(qX_long',mu_ini,sigma_ini.^2,w_ini),'--','color',colors.dre,'LineWidth',1),hold on
        title(['before ' Param.Optimizer{1,ixOpt} '; energy: ' num2str(Param.energy) ' MeV']); xlabel('depth')
    end
    
    switch currOpt
        case {'IPOPT'}
            addpath([pwd filesep 'optimization']);
            matRad_ipoptOptions
            options.lb                = vLowerBound;
            options.ub                = vUpperBound;
            
            options.ipopt.print_level = 4;
            options.ipopt.max_iter    = 15000;
            options.ipopt.tol         = 1e-13;
            
            funcs.objective           = @(x) matRad_APMfitObjective(x,qX_long,qY_long);
            funcs.gradient            = @(x) matRad_APMfitGradient(x,qX_long,qY_long);
            
            [result, ~]                = ipopt(wInit,funcs,options);
            
        case {'minimize'}
            
            returnRowVector = false;
            % minimze function by Carl Edward Rasmussen
            [result,~,~] = minimize(wInit', @(a) matRad_APMfitGausObjFunc(a,qX_long,qY_long,returnRowVector),2500);
            
        case {'fmincon'}
            options = optimset('Display', 'off','Algorithm','sqp','MaxIter',3000,...
                'MaxFunEvals',1e10,'TolFun',1e-10,'TolX',1e-9,'TolCon',1e-9,'GradObj','on');
            returnRowVector = true;
            result = fmincon(@(a) matRad_APMfitGausObjFunc(a,qX_long,qY_long,returnRowVector),wInit,[],[],[],[],vLowerBound,vUpperBound,[],options);
            
        case {'fminsearch'}
            
            options = optimset('Display', 'iter','FunValCheck','on','Algorithm','sqp','MaxIter',2000,...
                'MaxFunEvals',1e14,'TolFun',1e-14,'TolX',1e-25,'TolCon',1e-14,'GradObj','on');
            result = fminsearch(@(a) matRad_APMfitGausObjFunc(a,qX,qY,true),wInit,options);
            
    end
    
    [w_fit, mu_fit, sigma_fit] = splitOptResult(result);
    
    % evalute newly generate Gauss components
    mGaussComp = zeros(length(qX_long),Param.NumComp);
    for i = 1:Param.NumComp
        mGaussComp(:,i) = SG(qX_long,w_fit(i),mu_fit(i),sigma_fit(i));
    end
    
    % get newly fitted depth dose profile
    qYfitted = sum(mGaussComp,2)';
    
    % calcualte error
    if ixOpt == NumOpti
        Error = getError(qX_long,qY_long,qYfitted,Param);
    end
    
    if visBool
        Error = getError(qX_long,qY_long,qYfitted,Param);
        subplot(NumOpti,2,(ixOpt*2)),plot(qX_long,mGaussComp,'color',colors.dre),hold on;
        subplot(NumOpti,2,(ixOpt*2)),plot(qX_long,qY_long,'color',colors.red,'LineWidth',2);
        subplot(NumOpti,2,(ixOpt*2)),plot(qX_long,qYfitted,'--','color',colors.dkfzdB,'LineWidth',2), grid on, grid minor,
        title(['after ' Param.Optimizer{1,ixOpt} '; energy: ' num2str(Param.energy) ' MeV; max rel diff: ' num2str(Error.rel.maxDiff) ' %']);
        xlabel('depth'),set(gca,'Xlim',[min(qX) max(qX)]);
    end
    
    w_ini     = w_fit';
    mu_ini    = mu_fit';
    sigma_ini = sigma_fit';
    wInit     = [w_fit, mu_fit, sigma_fit];
    
end


% sanity check
vLowerBound(vLowerBound==1) = 0;
if sum(wInit<vLowerBound) > 0 || sum(wInit>vUpperBound) > 0
    warning('violations of boundaries'),
    ParamRecursion           = Param;
    ParamRecursion.Optimizer = {'IPOPT','fmincon'};
    [w_fit, mu_fit, sigma_fit, Error] = matRad_APMfitGausComp(qX,qY,ParamRecursion,[],true);
    
    w_ini     = w_fit';
    mu_ini    = mu_fit';
    sigma_ini = sigma_fit';
    wInit     = [w_fit, mu_fit, sigma_fit];
end



end




function Error = getError(x,yref,yhat,Param)

if isrow(yhat)
    yhat = yhat';
end

try
    % max squared error
    if isrow(yref)
        yref = yref';
    end
    diff = yref-yhat;
    Error.abs.maxSqError = max(diff.^2);
    % max deviation
    Error.abs.maxDiff     = max(abs(diff));
    
    % One way to define the relative difference of two numbers is to take their
    % absolute difference divided by the maximum absolute value of the two numbers.
    [val,ix]             = max(yref);
    rel_x                = x/x(ix);
    rel_y                = ((diff)/max([abs(yref)' abs(yhat)'])) * 100;
    Error.rel.vX         = linspace(0,1.5,600);
    Error.rel.vDiff      = interp1(rel_x,rel_y, Error.rel.vX,'linear',0);
    Error.rel.maxDiff    = max(abs(diff)/max([abs(yref)' abs(yhat)'])) * 100;
    Error.rel.meanDiff   = mean(abs(diff)/max([abs(yref)' abs(yhat)'])) * 100;
catch
    Error.rel.vX         = linspace(0,1.5,600);
    Error.rel.vDiff      = NaN;
    Error.rel.maxDiff    = NaN;
    Error.rel.meanDiff   = NaN;
    Error.abs.maxSqError = NaN;
    Error.abs.maxDiff   = NaN;
end

if Error.rel.maxDiff > 40
    warning(['rel diff bigger than 40, energyIx: ' num2str(Param.ixEnergy) '; NumCom: ' num2str(Param.NumComp)]);
end

end



function F = matRad_APMfitObjective(x,qX_long,qY_long)

[w, mu, sigma] = splitOptResult(x);

SG =  @(qX,qW,qMu,qSigma)((qW/(sqrt(2*pi*qSigma^2))).*exp(-((qX-qMu).^2)./(2*qSigma^2)));

qY = zeros(size(qY_long));

for i = 1:numel(w)
    qYini = SG(qX_long,w(i),mu(i),sigma(i));
    qY = qY + qYini;
end

F = sum((qY - qY_long).^2);

end

function g = matRad_APMfitGradient(x,qX_long,qY_long)
[ ~, g, ~ ] = matRad_APMfitGausObjFunc(x,qX_long,qY_long,true);
end



function [w, mu, ell] = splitOptResult(result)

% split the result into w mu ell
n = length(result);
if size(result,2) > 1
    w   = result(1:n/3);
    mu  = result(n/3+1:2*n/3);
    ell = result(2*n/3+1:end);
else
    w   = result(1:n/3)';
    mu  = result(n/3+1:2*n/3)';
    ell = result(2*n/3+1:end)';
end

if sum(w<-1e-3)>0
    warning('negative weights');
end
end
