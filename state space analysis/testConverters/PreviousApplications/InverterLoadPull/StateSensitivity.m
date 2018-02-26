function [ dXs ] = StateSensitivity( As, Bs, ts, u, varToPerturb, pI, dX, cI, HCM)
%   As, Bs, ts, u = state space matrix representation
%   varToPerturb = string for variable to perturb
%   pI = index of value to perturb in variable
%   dX = magnitude of perturbation
%   cI = index of variable to compensate (e.g. time compensation if
%       interval is perturbed to keep Ts)
%   HCM = half-cycle identity matrix (if used)

if(nargin == 7)
    cI = pI;
    dX2 = 0;
% elseif(nargin == 8)
%     dX2 = -dX;
end
dX2 = -dX;


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
    display('Variable to perturb not found');
    error('Variable to perturb not found');
end

if nargin < 9
    [ dXs] = SS_Soln( As, Bs, ts, u);
else
    [ dXs] = SS_Soln( As, Bs, ts, u, HCM);
end
    

end

