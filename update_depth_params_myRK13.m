%------------------------------------%
% UPDATE VALUES THAT DEPTH-DEPENDENT %
%------------------------------------%
function [iAttach, TempFunction]=update_depth_params_myRK13(z,r,sherwood_nu)
% Carbon dynamics Model per PARTICLE base

% Created by Trang Nguyen
% Modified: Apr 26, 2021

global varTemp TempCoeffArr TempAeArr TemprefArr Tkel 
global depth motility newAttachedCells

% Calculate TempFunction (temperature punishment)
    if varTemp
        Temp=12*exp(-z/150)+12*exp(-z/500)+2;
        Temp_100=12*exp(-100/150)+12*exp(-100/500)+2;%Temp at 100 m (original depth "depth")
        TempFun_100=TempCoeffArr.*exp(TempAeArr.*(1./(Temp_100+Tkel)-1./TemprefArr));
        TempFunction = TempCoeffArr.*exp(TempAeArr.*(1./(Temp+Tkel)-1./TemprefArr))/TempFun_100; %to make sure tempfun is 1 at 100m
    else
        TempFunction=1;
    end
        
 % Calculate attachment rate
    free_cell_function=exp(-3e-3*z)/exp(-3e-3*depth);
    iAttach(2)=4*pi*motility*r*sherwood_nu/100*free_cell_function*newAttachedCells*60*24; %cells/mins*24*60=cells/day, added to bopa2 pool
    iAttach(1)=0; %no new cells for original ba
    