function [ dXs ] = StateSensitivity(obj, varToPerturb, pI, dX, cI)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

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

[ dXs] = obj.SS_Soln();

    

end

