function stop = apm_fminconOutputFcn(x,optimValues,state)
stop = false;
switch state
    case 'iter'
        disp(['Search Direction: ' num2str(optimValues.fval)]);
end
end

