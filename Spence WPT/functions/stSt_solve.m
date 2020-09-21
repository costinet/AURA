function [Xs] = stSt_solve(As,Bs,ts,u)

% Commented by Spencer Cochran, October 22, 2019.

% -------------- Steady State Solver for Spencer --------------
% ----------------------- Jared Baxter ------------------------
% --------------------- October 14, 2019 ----------------------

n = size(As,3);     % # of intervals. 
ns = size(As,1);    % # of states. 
pi = length(u);     % # of inputs. 
    
for i = 1:1:n                                               % For each interval.
    Aa(:,:,i) = [As(:,:,i), Bs(:,:,i); zeros(pi,ns+pi)];    % Augmented: Aa = [A, B; zeros]; 
    expAa(:,:,i) = expm(Aa(:,:,i)*ts(i));                   % e^(Aa*t)
        if (i==1) expAa_mult = expAa(:,:,1);                % Multiply all intervals. 
        else      expAa_mult = expAa(:,:,i)*expAa_mult; end
end

x0b = null(eye(ns+pi) - expAa_mult);	% Basis for steady state solution. 
Ks = x0b(ns+1:end, :) \ u;            	% Matrix Left Division: calculate basis scaling factor.
Xs = x0b*Ks;                         	% Scale basis for final augmented steady state. 

for i = 1:1:n
    Xs(:,end+1) = expAa(:,:,i)*Xs(:,i); end

assert(norm(Xs(:,1)-Xs(:,end)) < 1e-9); %Check that steady-state was solved correctly
% -------------------------------------------------------------

Xs(end-pi+1:end,:) = [];    % Spencer Adds: Remove the input vector from augmentation.

end