function color = colorspecs()

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

color.lightdre = [1,1,1] - 0.5 * ([1,1,1] - color.dre);
color.lightblu = [1,1,1] - 0.5 * ([1,1,1] - color.blu);
color.lightora = [1,1,1] - 0.5 * ([1,1,1] - color.ora);
color.lightgre = [98/255 202/255 4/255];

color.darkgre = [37/255 111/255 19/255];
color.darkora = [196/255 137/255 0/255];
color.orared = [218/255 89/255 33/255];

end