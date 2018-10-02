function covDose = apm_calcCovDose(Vijlm,w)

covDose = reshape(reshape(permute(Vijlm,[1 3 4 2]),size(Vijlm,1)^2,numel(w)^2) * reshape(w*w',numel(w)^2,[]),size(Vijlm,1),size(Vijlm,1)); %Manual product via matricication
%covDose =  double(ttt(tensor(Vijlm),tensor(w*w'),[4 2],[1 2])); %Product with tensor toolbox

end

