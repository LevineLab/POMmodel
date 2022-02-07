%% Carbon dynamics Model per PARTICLE base %% 
% Model input parameters % 

% Created by Trang Nguyen
% Modified: Apr 26, 2021

%% Global variables
global nbft initba c_avg_pa rcell
global mlin D kmom ymom vmom_max varTemp
global depth L0 newAttachedCells cell_c_content cell_per_particle
global day_to_sec motility
global Af D_Omand D_s c_avg_pa_ref fluid_density_ref g mu rs_3_minus_d cut_off volu convert_Alldredge_Omand %For Omand density and sinking equation
global viscosity sc
global betas Ws1
global  ndays  timestep dt TempFun
global TempAeArr TemprefArr Tkel TempCoeffArr
%-------------------------------------------------------%
%%%% LOAD IN USER DEFINED PARAMETERS
%-------------------------------------------------------%

Input_Parameters

%-------------------------------------------------------%
%           MODEL PARAMS                               %
%-------------------------------------------------------%

% Constants
day_to_sec=24*60*60;
D=6.7E-6*3600*24;%cm2.day-1 %Diffusion of monomer (D of glucose)
viscosity=1E-2;%cm2.s-1 %Kioebe 2001 
sc=viscosity/(D/3600/24); %Schmidt number

%Depth and Particle Radius
depth=100; %m  Depth of POM formation
rcell=6.2e-5; %cm  for Vol_cell=1um3

%% POM density and sinking rate calculation %%
Af=1.49; %unitless
D_Omand=1.85;%1.85; %aggregate; =3 for small solid particle
D_s=3; %small particle
c_avg_pa_ref=1.23; %g/cm3, McDonell 2010, Sicko-Goad 1980, Cram 2018 - Measured pom density for small solid particles
fluid_density_ref=1.028;%g/cm3, Omand 2020, at the base of euphotic zone, assumed to be constant in the density equation
g=9.81; %m/s2
mu=1e-3; %kg/m/s
rs_3_minus_d=0.0012*18*mu/(g*Af*(c_avg_pa_ref-fluid_density_ref)*1000*(1e-3)^(D_Omand-1));%rs is in meter; in which the reference particle is 1e-3 m with 0.0012 m/s sinking rate
cut_off=rs_3_minus_d^(1/(3-D_Omand))*100;%cm, aka radius of the sub-unit
volu=4/3*pi*r^3;

f=Af*rs_3_minus_d*(r/100).^(D_Omand-3);%unitless, Omand 2020, Fraction of POM in an aggregate
[I]=find(f>1);
f(I)=1;
[I]=find(r<=cut_off);
f(I)=1;

c_avg_pa_omand=f*(c_avg_pa_ref-fluid_density_ref)+fluid_density_ref; %g/cm3, Omand 2020, POM Note4 Trang
[I]=find(r<=cut_off);
c_avg_pa_omand(I)=c_avg_pa_ref;%Assuming all r<100um particles have constant density
[I]=find(c_avg_pa_omand>=c_avg_pa_ref);
c_avg_pa_omand(I)=c_avg_pa_ref;%Assuming all particles > constant density to have constant density
pseudo_pom=f.*c_avg_pa_ref.*volu/12*1000;%mmolC/particle, Omand 2020, POM Note4 Trang
convert_Alldredge_Omand=0.0494;%hard code the B ~ regression between Omand C content and Alldredge C content to account for non-org C fraction in a particle
c_avg_pa=pseudo_pom*convert_Alldredge_Omand;%save this B: 0.0494

w_aggregate=g*Af*(c_avg_pa_ref-fluid_density_ref)*1000*(r/100).^(D_Omand-1)*rs_3_minus_d/(18*mu)*100;%cm/s, Omand 2020, POM Note4 Trang, Sinking rate based on density at certain depth and based on the change of fluid density at depth
w_smallpa=g*Af*(c_avg_pa_ref-fluid_density_ref)*1000*(r/100).^(D_s-1)/(18*mu)*100;%cm/s, Sinking rate for solid particle fraction dimension D=3 aka a_s^(D-3)=1
if r<=cut_off
    Ws_cms=w_smallpa; %cm/s
else
    Ws_cms=w_aggregate;%cm/s
end
Ws=Ws_cms*864; %m/day
Ws1=Ws;

rey= Ws_cms * r / viscosity; %unitless, Reynold number %Reynold Number (Moradi 2018, describe how viscous the particle is compared to the environment, depends on the diffusion and advection on the particle %Re=sinking*radius/water viscosity (high Re=High sinking)
sc=viscosity/(D/3600/24); %Schmidt number, unitless, Bianche 2018 Nat Geo (Supp)
sh=1+0.619*rey^0.412*sc^(1/3);%Sherwood number, unitless, Bianche 2018 Nat Geo (Supp)
diffloss=4*pi*D*sh*r*1e-6;%m3/day/particle

%% Bacterial initialization

%Initial population
nbft=2;%number of bacterial groups so that initial colonizers and 'recruits' tracked seperately
cell_c_content=20E-15/12*1E3; %mmolC.cell-1 %Lee and Fuhrman 1987, unit fg into mmol C
initba=cell_per_particle*cell_c_content;%mmolC.particle-1

%Encounter rate
newAttachedCells=F0_max*cell_c_content;%mmolC of cell/m3


%% Bacterial dynamics and behavior - Global variables

% Motility (For Encounter Rate calculation)
motility=(1e-9 + 1e-10)/2*60; %m2/minute

%Mortality rates
mlin=2.5;%linear mortality for b

%Detachment (Leaving) rate
L0=0.25;

%Metabolism for Aerobic Heterotrophs
ymom=0.2; %mmolC_cells/mmolC_MOM - aerobic bacteria yield
vmom_max=umax/ymom;%day-1
kmom=7.2;%mmolC_MOM.m-3; Half saturation for MOM

%% Initialize following recruitment experiment
betas=repmat( mybeta, [1,nbft] );
bo_pas(1)=initba;%early colonizers
bo_pas(2)=cell_c_content; %recruited cells ~ late comers

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%
%-------------------------------------------------------%
%           PARAMS REGARDING TO SOLVING EQUATIONS       %
%-------------------------------------------------------%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%
%% TimeSteps
ndays=5;%day
dt=1/24/6;%1/24/6/1; %day/step ~ apx every 10 mins
timestep=round(ndays/dt);

%% Temperature Function through depth
varTemp=1;   %1 means variable temp, 0 means constant temp
TempAeArr = -4E3;  
TemprefArr = 293.15 ;      
Tkel = 273.15 ;  
TempCoeffArr = 0.8;

Temp=12*exp(-depth/150)+12*exp(-depth/500)+2;
TempFun = TempCoeffArr.*exp(TempAeArr.*(1./(Temp(1)+Tkel)-1./TemprefArr));

if varTemp==0
    TempFun=1;
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%
%-------------------------------------------------------%
%           INITIALIZE VALUES FOR PATCH                %
%-------------------------------------------------------%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%

PATCH.pom=c_avg_pa; %mmol.particle-1
PATCH.mom=0; %mmol.m-3.particle-1
PATCH.z=depth;
%PATCH.pomz=0;
for ibft=1:nbft
    PATCH.bo_pa(ibft,1)=bo_pas(ibft);
    PATCH.celldensity(ibft,1)=cell_density_init/nbft;%cells.m-2.particle-1
    PATCH.beta(ibft,1)=betas(ibft);%day-1
end

PATCH.pool=nbft;
PATCH.depth=depth;
PATCH.time=0;
PATCH.r=r;
% PATCH.TemF=TempFun;
PATCH.newcells=newAttachedCells;
PATCH.varTemperature=varTemp;
PATCH.vmax=vmom_max*ymom-mlin;
PATCH.km=kmom;
PATCH.mlinear=mlin;
