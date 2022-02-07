%Input Parameters

r=0.05;%cm 

mybeta=150;

%you can either set total cells per particle or cell density at formation
cell_per_particle=1.5E4; %initial cells particle-1 
cell_density_init=cell_per_particle/(4*pi*r*r*1E-4); %cells.m-2.particle-1

F0_max=5e10;%5e10; %cells/m3  concentration of free living bact community 

umax= 3.7;  %day-1

runName=['test'];