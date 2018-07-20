function [w_ini, mu_ini, sigma_ini] = matRad_APMfitGetIniParameters(qX,qY,Param)

% this function defines manually iniital fitting parameter

SG =  @(qX,qW,qMu,qSigma)((qW/(sqrt(2*pi*qSigma^2))).*exp(-((qX-qMu).^2)./(2*qSigma^2)));

if strcmp(Param.Type,'LET')
    ix  = find(qY >max(qY)*0.8,1,'first');
    if ix == numel(qX)
        ix = ix -1;
    end
    val = qY(ix);
else
    [val,ix] = max(qY);
end

peakPos  = qX(ix);

switch Param.radMod
    
    case {'protons'}
        
        mu_ini    = peakPos * linspace(0.05,0.9,Param.NumComp)';
        sigma_ini = 1*(peakPos/Param.NumComp)* ones(Param.NumComp,1);
        
        val = qY(1)/SG(0,1,0,sigma_ini(1));
        if strcmp(Param.Type,'LET')
            w_ini = interp1([qX(1) qX(ix)*0.5 qX(ix)*0.75 qX(ix) qX(end)],[val*0.4 val*0.7 val*3 val*10 val*12],mu_ini) *  1 +...
                (Param.vEnergies(Param.ixEnergy)/Param.vEnergies(end))*60;
        else
            w_ini = interp1([qX(1) qX(ix)*0.5 qX(ix)*0.75 qX(ix) qX(end)],[val*0.8 val*1 val*1.3 val*3 val*0.8],mu_ini) *  1 +...
                (Param.vEnergies(Param.ixEnergy)/Param.vEnergies(end))*60;
        end
        
    case {'carbon'}
        
        if isequal(Param.Type,'Z')
           
            mu_ini    = peakPos * linspace(0.02,1.05,Param.NumComp)';
            sigma_ini = 0.5*(peakPos/Param.NumComp)* ones(Param.NumComp,1);
            val = qY(1)/SG(0,1,0,sigma_ini(1));
            w_ini = interp1([qX(1) qX(ix)*0.5 qX(ix)*0.75 qX(ix) qX(end)],[val*0.8 val*0.85 val*1 val*3 val*0.8],mu_ini);% *...
            
        elseif isequal(Param.Type,'alphaDose')
            
            mu_ini    = peakPos * linspace(0.05,1.15,Param.NumComp)';
            Fac       = peakPos/30;
            sigma_ini = interp1([qX(1) qX(ix) qX(end)],[Fac Fac*0.75 Fac*100],mu_ini);
            
            try
                w_ini = interp1([qX(1) qX(ix) qX(end)],[Fac*val*0.25 Fac*val*0.75 Fac*val*5],mu_ini);
            catch
                error(['Weight INI error at ' num2str(peakPos)])
                w_ini = 30*ones(Param.NumComp,1);
            end
            
        elseif isequal(Param.Type,'SqrtBetaDose')
            
            mu_ini    = peakPos * linspace(0.05,1.11,Param.NumComp)';
            Fac       = peakPos/30;
            sigma_ini = interp1([qX(1) qX(ix) qX(end)],[Fac Fac*0.75 Fac*100],mu_ini);
            try
                w_ini = interp1([qX(1) qX(ix) qX(end)],[Fac*val*1.25 Fac*val*1.75 Fac*val*5],mu_ini);
            catch
                error(['Weight INI error at ' num2str(peakPos)])
                w_ini = 30*ones(Param.NumComp,1);
            end
            
        elseif isequal(Param.Type,'LET')
            
            Param.NumComp = Param.NumComp-2;
            
            mu_ini    = peakPos * linspace(0.05,1.1,Param.NumComp)';
            sigma_ini = 0.8*(peakPos/Param.NumComp)* ones(Param.NumComp,1);
            
            val = qY(1)/SG(0,1,0,sigma_ini(1));
            if strcmp(Param.Type,'LET')
                w_ini = interp1([qX(1) qX(ix)*0.5 qX(ix)*0.75 qX(ix) qX(end)],[val*0.4 val*0.7 val*3 val*10 val*12],mu_ini) *  1 +...
                    (Param.vEnergies(Param.ixEnergy)/Param.vEnergies(end))*60;
            else
                w_ini = interp1([qX(1) qX(ix)*0.5 qX(ix)*0.75 qX(ix) qX(end)],[val*0.8 val*1 val*1.3 val*3 val*0.8],mu_ini) *  1 +...
                    (Param.vEnergies(Param.ixEnergy)/Param.vEnergies(end))*60;
            end
            
            mu_ini(end+1)    = peakPos * 1.3;
            sigma_ini(end+1) =  sigma_ini(end)*1.5;
            w_ini(end+1)     = w_ini(end)/2;
            
            mu_ini(end+1)    = peakPos * 1;
            sigma_ini(end+1) =  sigma_ini(end)*20;
            w_ini(end+1)     = w_ini(end);
            
            Param.NumComp = numel(w_ini);
        end
        
        
end



end

