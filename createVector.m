%--------------------------------------------------%
% Modified createVector for Trang simplified model %
%                Not working yet                   %
%--------------------------------------------------%

function [dCdt]=createVector(NPATCH)
global nbft

dCdt(1)= NPATCH.pom(end); %total pool
dCdt(2)= NPATCH.mom(end); %total pool
dCdt(3)= NPATCH.z(end); %depth over time inside myRK
%dCdt(4)= NPATCH.pomz(end); %POM over depth over time inside myRK

for ibft=1:nbft 
    dCdt(3+ibft)= NPATCH.bo_pa(ibft,end);%dCdt(4+ibft)= NPATCH.bo_pa(ibft,end);
end
