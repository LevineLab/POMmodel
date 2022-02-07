%--------------------------------%
% UPDATE VALUES THAT R-DEPENDENT %
%--------------------------------%
function [Ws_rk,diffloss,sherwood_nu] = update_r_values_myRK13(r)

% Carbon dynamics Model per PARTICLE base

% Created by Trang Nguyen
% Modified: Apr 26, 2021

global viscosity sc D 
global Af D_Omand D_s c_avg_pa_ref fluid_density_ref g mu rs_3_minus_d cut_off  %For Omand density and sinking equation

w_aggregate=g*Af*(c_avg_pa_ref-fluid_density_ref)*1000*(r/100).^(D_Omand-1)*rs_3_minus_d/(18*mu)*100;%cm/s, Omand 2020, POM Note4 Trang, Sinking rate based on density at certain depth and based on the change of fluid density at depth
w_smallpa=g*Af*(c_avg_pa_ref-fluid_density_ref)*1000*(r/100).^(D_s-1)/(18*mu)*100;%cm/s, Sinking rate for solid particle fraction dimension D=3 aka a_s^(D-3)=1
if r<=cut_off
    Ws_cms=w_smallpa; %cm/s
else
    Ws_cms=w_aggregate;%cm/s
end
Ws_rk=Ws_cms*864; %m/day

rey=Ws_cms*r/viscosity; %unitless
sherwood_nu=1+0.619*rey^0.412*sc^(1/3);%Sherwood number, unitless, Bianche 2018 Nat Geo (Supp)

diffloss=4*pi*D*sherwood_nu*r*1e-6;%m3/day/particle
