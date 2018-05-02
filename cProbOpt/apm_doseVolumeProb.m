function [expDV,stdDV,stdDV_ub] = apm_doseVolumeProb(expDose,covDose,dParam,method)

if nargin < 4
    method = 'int_gauss';
end


switch method
    case 'int_gauss'
        %summands = 1 - normcdf(dParam* ones(numel(expDose),1),expDose,sqrt(diag(covDose)));
        %exp = 1/numel(expDose) * sum(summands);
        %expDV = 1 / (2*numel(expDose)) * sum(1 -  erf((dParam-expDose)./sqrt(2*diag(covDose))));
        %exp = 1 - mvncdf(dParam*ones(numel(expDose),1),expDose,covDose);
        %exp = numel(expDose)/2 * sum(erfEval);
        %Fval = normcdf(dParam * ones(numel(expDose),1),expDose,sqrt(diag(covDose)));
        Fval = 0.5*(1 + erf((dParam * ones(numel(expDose),1) - expDose)./sqrt(2*diag(covDose))));
        %Fval = heaviside((dParam * ones(numel(expDose),1) - expDose) ./ sqrt(diag(covDose)));
        expSum = sum(1 - Fval);
        expDV = 1/numel(expDose) * expSum;
        stdDV_ub = sqrt(1/numel(expDose) .* sum((1-Fval).*Fval));
        
        biCum = 0;
        
        for i = 1:numel(expDose)
            sig_i = sqrt(covDose(i,i));
            for l = i+1:numel(expDose)
                sig_l = sqrt(covDose(l,l));
                r = covDose(i,l) / (sig_i*sig_l);
                prob = bvn((dParam - expDose(i)) / sig_i,Inf,(dParam - expDose(l)) / sig_l,Inf,r);
                %prob = mexBVNcdf([expDose(i),expDose(l)],[dParam,dParam],[covDose(i,i),covDose(i,l); covDose(l,i) covDose(l,l)]);
                biCum = biCum + prob;
                %{
                xi = linspace(expDose(i) -  3*sig_i,expDose(i) + 3*sig_i,50);
                xl = linspace(expDose(l) -  3*sig_l,expDose(l) + 3*sig_l,50);
                [Xi,Xl] = meshgrid(xi,xl);
                F = mvnpdf([Xi(:) Xl(:)],expDose([i l])',covDose([i l],[i l]));
                F = reshape(F,length(xi),length(xl));
                
                
                imagesc(xi,xl,F); colormap(flipud(gray)); colorbar;
                hold on;
                set(gca,'YDir','normal');
                yli = get(gca,'YLim');
                xli = get(gca,'XLim');
                plot([dParam dParam],[yli(1) yli(2)],'r');
                plot([xli(1) xli(2)],[dParam dParam],'r');
                title(['Probability d > dT: ' num2str(prob) ]);
                waitforbuttonpress;
                hold off;
                %}
            end
        end
                
        stdDV = 1/numel(expDose)^2 * (expSum + 2*biCum);
        
        %stdDV = biCum / numel(expDose)^2;
        
        stdDV = sqrt(stdDV - expDV^2);
        
    otherwise
        error(['Method ''' method ''' not defined!']);
end

end

