function [dXs] = Baxter_StateSensitivity(obj, varToPerturb, pI, dX, cI)
%Baxter_StateSensitivity Solves the state space for a slightly perturbed variable and solve for how that affects the steady state of the converter
%   varToPerturb is on of the char arrays 'As','Bs','ts','u'
%
%   Pi sets the increase in dX value that is desired, if As or Bs then
%   it needs to be a 1x3 double of the position of the value, if ts or
%   u then a 1x1 double to indicate the time period or the contestant
%   value index
%
%   dX is the amount by which the value referenced by varToPerturband
%   and pI should be perturbed
%
%   cI is the variable that should be reduced in order to maintain
%   constant energy in the circuit or switching frequency
%
%   For example: if time interval 2 is perturbed then time must be
%   taken out of some other time interval in the circuit or a change
%   in switching frequency will result

%     _   _   _  ____    _    
%    / \ | | | |/ _  |  / \   
%   / _ \| | | | (_| | / _ \  
%  / ___ | |_| |> _  |/ ___ \ 
% /_/   \_\___//_/ |_/_/   \_\


As = obj.As;
Bs = obj.Bs;
ts = obj.ts;
u = obj.u;

if(nargin == 4)
    cI = pI;
    dX2 = 0;
elseif(nargin == 5)
    dX2 = -dX;
end

if(varToPerturb(1:2) == 'As')
    As(pI(1), pI(2), pI(3)) = As(pI(1), pI(2), pI(3)) + dX;
    As(cI(1), cI(2), cI(3)) = As(cI(1), cI(2), cI(3)) + dX2;
elseif(varToPerturb(1:2) == 'Bs')
    Bs(pI(1), pI(2), pI(3)) = Bs(pI(1), pI(2), pI(3)) + dX;
    Bs(cI(1), cI(2), cI(3)) = Bs(cI(1), cI(2), cI(3)) + dX2;
elseif(varToPerturb(1:2) == 'ts')
    ts(pI) = ts(pI) + dX;
    ts(cI) = ts(cI) + dX2;
elseif(varToPerturb(1) == 'u')
    u(pI) = u(pI) + dX;
    u(cI) = u(cI) + dX2;
else
    disp('Variable to perturb not found');
    error('Variable to perturb not found');
end

[dXs] = obj.SS_Soln(0,As,Bs,ts,u); % Zero is for keeping Xss from last time, not sure if this is applicable here. might want to pass a variable in the future
[dXs] = obj.CorrectXs(0,dXs);


end
