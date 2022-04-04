function c = apm_calcCopula(u1, u2, par, copula_model)
% Returns Copula(u1,u2; par), using the copula model in the input,
% where 'Copula' is the joint cummulative distribution function of u1 and u2
% with parameter par.

switch true
    case (strcmpi(copula_model,'gauss') == 1) % Gaussian Copula
        r = par; % Correlation coefficient
        v1 = norminv(u1); % Inverse Gaussian CDF
        v2 = norminv(u2);
        c = arrayfun(@(v1,v2,r) bvn(-Inf,v1,-Inf,v2,r), v1, v2, r); 
    
    case (strcmpi(copula_model,'amh') == 1) % Ali-Mikhail-Haq
        c = u1 .* u2 ./ (1 - par .* (1-u1) .* (1-u2));
    
    case (strcmpi(copula_model, 'gumbel') == 1) % Gumbel
        c = copulacdf('Gumbel', [u1(:),u2(:)], par);
        c = reshape(c,numel(u1(:,1)), numel(u2(:,2)));
    
    case (strcmpi(copula_model, 'clayton') == 1) % Clayton
        c = copulacdf('Clayton', [u1(:),u2(:)], par);
        c = reshape(c,numel(u1(:,1)), numel(u2(:,2)));
    
    case (strcmpi(copula_model, 'frank') == 1) % Frank
        c = copulacdf('frank', [u1(:),u2(:)], par);
        c = reshape(c,numel(u1(:,1)), numel(u2(:,2)));
    %case (strcmpi(copula_model, 'uniform') == 1) % Uniform perfect correlation
    %    c = min(u1, u2);
    %    if par ~= 1
  	%        disp(['Unapropriate Copula model!']);
    %    end
    otherwise
        disp (['Invalid Copula model']);
end

end
