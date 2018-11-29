function [t] = deadtimecalc(obj,start_value,end_value,state_pos,time_pos)
%deadtimecalc calculates the deadtime needed for a given A and B
%matrix and an associated time inverval


As = obj.As;
Bs = obj.Bs;
Xs = obj.Xs;
u = obj.u;
ts = obj.ts;

%% Slightly Decrease diagonal elements to create invertable matrix
% (Within MATLAB bounds)
% Similar principle to SS_Soln.m for finding the natural response of
% the circuit.

hi = -6; % Intial guess (might be able to educated guess this)
ns = size(As,1);
if cond(As(:,:,time_pos))>1*10^9 % 10^9 is an educated guess
    invA = (As(:,:,time_pos)-(1-10^(hi))*eye(ns))^-1; % Calculate the educated guess value
    
    %%%% Use cond() to see where how 'able' the matrix is to converge %%%%
    %%% Can also try and use optimization commented out below %%%%
    
    % Loop through increaseing hi until the condition of the inverted
    % matrix reaches a 'large' value
    while cond(As(:,:,time_pos)-(1-10^(hi))*eye(ns))<1*10^9
        hi = hi-1;
        invA = (As(:,:,time_pos)-(1-10^(hi))*eye(ns))^-1; % Calculate the educated guess value
    end
else
    invA = (As(:,:,time_pos)-(1-10^(hi))*eye(ns))^-1; % Calculate the educated guess value
end

invA_select = invA(state_pos,:);
B = Bs(:,:,time_pos);
t = (log(end_value + invA_select*B*u)-log(start_value + invA_select*B*u))*invA_select;



end

