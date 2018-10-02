function ttAvailable = apm_tensorToolboxAvailable()
    ttAvailable = exist('tensor','class') && exist('ttv','function') && exist('ttm','function'); 
end

