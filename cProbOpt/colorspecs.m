function color = colorspecs()

color.gre1 = [0,0.4717,0.4604]; 
color.dre = [0.4906,0,0]; 
color.ora = [255,153,51] ./ 255;
color.blu = [0,0,0.509];
color.gra = 0.5 * ones(3,1);
color.graDark = 0.2 * ones(3,1);
color.red = [188/255,18/255,61/255];
color.gre = [1/255,126/255,51/255];

color.dkfzlB = [169/255,188/255,213/255]; 
color.dkfzmB = [70/255,118/255,173/255]; 
color.dkfzdB = [0, 75/255, 142/255]; 


color.lightmpg = [1,1,1] - 0.5 * ([1,1,1] - color.gre1);
color.lightdre = [1,1,1] - 0.5 * ([1,1,1] - color.dre);
color.lightblu = [1,1,1] - 0.5 * ([1,1,1] - color.blu);
color.lightora = [1,1,1] - 0.5 * ([1,1,1] - color.ora);
color.lightgre = [98/255 202/255 4/255];
color.lightgre = [98/255 202/255 4/255];
color.lightgre2 = [198/255 250/255 143/255];
color.lightgre3 = [210/255 255/255 169/255];

color.darkgre = [37/255 111/255 19/255];
color.darkora = [196/255 137/255 0/255];
color.orared = [218/255 89/255 33/255];

color.plotblue = [0.4549,0.6784,0.8196];
color.plotorange = [0.8431,0.1882,0.1529];
color.plotyellow = [0.9922,0.6824,0.3804];

color.plotblue_light = [1,1,1] - 0.5 * ([1,1,1] - color.plotblue);
color.plotorange_light = [1,1,1] - 0.5 * ([1,1,1] - color.plotorange);
color.plotyellow_light = [1,1,1] - 0.5 * ([1,1,1] - color.plotyellow);

color.plotblue_grayed = hsv2rgb([1 0.25 1] .* rgb2hsv(color.plotblue));
color.plotorange_grayed = hsv2rgb([1 0.25 1] .* rgb2hsv(color.plotorange));
color.plotyellow_grayed = hsv2rgb([1 0.25 1] .* rgb2hsv(color.plotyellow));

color.plotblue2 = [0.2706,0.4588,0.7059];
color.plotorange2 = [0.9569,0.4275,0.2627];
color.plotyellow2 = [0.4549,0.6784,0.8196];

color.plotblue2_light = [1,1,1] - 0.5 * ([1,1,1] - color.plotblue2);
color.plotorange2_light = [1,1,1] - 0.5 * ([1,1,1] - color.plotorange2);
color.plotyellow2_light = [1,1,1] - 0.5 * ([1,1,1] - color.plotyellow2);

end