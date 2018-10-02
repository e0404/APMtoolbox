function sqDev = apm_squaredDev(dose,dParam)

sqDev = sum((dose-dParam).^2);

end

