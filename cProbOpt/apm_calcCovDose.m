function covDose = apm_calcCovDose(Vijlm,w)

if isa(Vijlm,'sptensor')
    %covDose = double(ttt(Vijlm,sptensor(w*w'),[4 2],[1 2]));
    covDose = ttv(Vijlm,w,4);
    covDose = ttv(covDose,w,2);
    covDose = double(covDose);
else
    covDose = reshape(reshape(permute(Vijlm,[1 3 4 2]),size(Vijlm,1)^2,numel(w)^2) * reshape(w*w',numel(w)^2,[]),size(Vijlm,1),size(Vijlm,1)); %Manual product via matricication
end

end

