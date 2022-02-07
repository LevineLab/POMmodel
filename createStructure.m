%-----------------------------------------------------%
% Modified createStructure for Trang simplified model %
%                Not working yet                      %
%-----------------------------------------------------%

function [NPATCH]=createStructure(NPATCH, dCdt, t)
global nbft

NPATCH.pom(t)=dCdt(1);%total pool
NPATCH.mom(t)=dCdt(2);%total pool
NPATCH.z(t)=dCdt(3);%depth over time, inside myRK
%NPATCH.pomz(t)=dCdt(4);%POM over depth over time, inside myRK

for ibft=1:nbft 
    NPATCH.bo_pa(ibft,t)=dCdt(3+ibft);%dCdt(4+ibft);
end
