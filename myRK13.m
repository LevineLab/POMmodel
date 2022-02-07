function [dCdt] = myRK13(t,C)
% Carbon dynamics Model per PARTICLE base

% Created by Trang Nguyen
% Modified: March 25, 2021

global vmom_max ymom kmom
global mlin L0 nbft betas cell_c_content
global dt Ws1
global Af D_Omand c_avg_pa_ref rs_3_minus_d cut_off convert_Alldredge_Omand %For Omand density and sinking equation

%% Initialize the Vector
pom=C(1); %mmolC.particle-1
mom=C(2); %mmolC.m-3
z=C(3); %m

for ibft=1:nbft 
    bo_pa(ibft)= C(3+ibft);%C(4+ibft); %mmolC.particle-1
end

%% Update values with the change of particle radius during sinking
new_c_avg_pa=pom;

f_mult_vol=new_c_avg_pa/convert_Alldredge_Omand/(c_avg_pa_ref/12*1000);
pseudo_r1=(3/4*f_mult_vol/pi)^(1/3); %r calculated if f=1 aka small particles
pseudo_r2=(f_mult_vol/(4/3*pi*Af*rs_3_minus_d*100^(3-D_Omand)))^(1/D_Omand); %r calcualted if f<1  
if pseudo_r1<= cut_off
    my_r=pseudo_r1;
else
    my_r=pseudo_r2;
end
  
[Ws_rk,diffloss,sherwood_nu] = update_r_values_myRK13(my_r);

%% Update params with the change of depth 

[iAttach, TempFunction]=update_depth_params_myRK13(z,my_r,sherwood_nu); %iAttach is cell/day
    
%% Bacterial Growth for each bacteria group
%POM dynamics
tot_ba=sum(bo_pa);
for ibft=1:nbft
    dpomsdt(ibft)= -betas(ibft) * bo_pa(ibft);
end
    
%MOM dynamics
tot_mom=mom;
tot_vmom=vmom_max.*tot_mom./(tot_mom+kmom).*TempFunction; %Rate of monomer uptake
tot_u_bo_pa=tot_vmom*ymom; %Bacteria growth rates
tot_bact_off=L0;%Bacterial rate of leaving particle

% Calculate pseudo steady state of MOM
% Solve for MOM when 0=Production-Consumption-Diffusion 
a=diffloss;%m3/day/particle %Diffusion
b=tot_ba.*vmom_max; %Consumption
p=tot_ba.*betas(1); %beta(1) and beta(2) are the same
if a==0
    x=kmom*p/(b-p);
else
   x=sqrt((kmom*p+a*((a*kmom+b-p)/(2*a))^2)/a)-(a*kmom+b-p)/(2*a); %Pseudo-Steady state: mom(i+1)
end
%End of pseudo-steady state calculation  

for ibft=1:nbft
   bact_off(ibft)=tot_bact_off;
   vmom(ibft)=tot_vmom;
   u_bo_pa(ibft)=vmom(ibft)*ymom;
end 

for ibft=1:nbft
   dbo_padt(ibft)= bo_pa(ibft).*(tot_u_bo_pa - mlin.*TempFunction - tot_bact_off);
end 
   dmomdt=(x-tot_mom)/dt; %Place holder, not used - just back calculated from Steady state approximation

%% Derivative equations

for ibft=1:nbft
         dCdt(3+ibft)= dbo_padt(ibft)+iAttach(ibft);
%        dCdt(4+ibft)= dbo_padt(ibft)+iAttach(ibft);
end

dpomdt=dpomsdt(1); %Ba group 1

if ibft>1
    for ibft=2:nbft
        dpomdt=dpomsdt(ibft)+dpomdt;
    end
end

dzdt=Ws_rk;% m/day

dCdt(1)= dpomdt;%total pool of POM
dCdt(2)= dmomdt;%total pool of MOM
dCdt(3)= dzdt;%depth

dCdt=dCdt';
