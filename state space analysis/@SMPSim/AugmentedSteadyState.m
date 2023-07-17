function [ Xs] = AugmentedSteadyState(obj, dts)
% Steady-state solution of switched system using state-space matrices
%
% [ Xs] = SS_Soln( As, Bs, ts, u, Xi, Bi) finds the state values Xs in
% steady-state for the switched system described by As, Bs, and ts
% according to
%
% dx/dt = Ai*x(t) + Bi*u ,
%
% for the ith interval.
%
% As is a 3-dimensional matrix of values for Ai, where As(:,:,i) is the 2D
% square matrix Ai during the ith interval.  
% Bs is a 3-dimensional matrix of values for Bi, where Bs(:,:,i) is the 2D
% matrix/vector Bi during the ith interval.
% ts is a vector of the time durations of each inverval
% u is the (assumed constant) independent input vector
% Xi is a vector "guess" of the results, which is used to determine
% solution validity in the case of singular matrices
% Bi is a vector of boundary coefficiencts used when iteration is required
% due to non-convergence of a solution.  A valid solution is required to be
% bounded within Bi.*Xi <= Xs <= (2-Bi).*Xi
%
% The result, Xs is a 2D matrix with ns rows and n+1 columns, where ns is the
% number of states in the system and n is the number of switching
% intervals.  The first and last column of Xs should be identical,
% corresponding to a valid steady-state solution.

    if nargin <2
        dts = obj.ts*0;
    end

    As = obj.As;
    Bs = obj.Bs;
    ts = obj.ts + dts;

    if isempty(obj.IHC)
        obj.IHC = eye(size(As,1));
    end
    
    if size(obj.u,3) == 1
        %% HEY THIS SHOULD BE MOVED ELSEWHERE
        u = repmat(obj.u,1,1,(length(ts)));
    elseif size(obj.u,3) == length(ts)
        u = obj.u;
    else
        error('invalid input vector u specificed');
    end

    n = size(As,3);
    ns = size(Bs,1);

    In = eye(size(As,1) + 1, size(As,2) + 1);
    EA = In;
    Atil = zeros(size(EA,1), size(EA,2), n);
    for i = n:-1:1
        Atil(:,:,i) = [As(:,:,i), Bs(:,:,i)*u(:,:,i); zeros(1, size(EA,1))];
        expAtil(:,:,i) = expm(Atil(:,:,i)*ts(i));
        EA = EA*expAtil(:,:,i);
    end

    IHC = eye(size(EA));
    IHC(1:ns,1:ns) = obj.IHC;


    depStates = sum(abs(eye(ns)-obj.Is(:,:,1)),2)  ~= 0;
    depStates(end+1) = 0;
    nd = sum(depStates);

    Zn = null(In(~depStates,~depStates)-IHC(~depStates,~depStates)*EA(~depStates,~depStates));

    k = Zn(end,:)\1;
    Xss = zeros(ns,1);
    Xss(~depStates(1:end-1)) = Zn(1:end-1,:)*k;
    Xss = obj.Is(:,:,1)*Xss;

%     EA*[Xss; 1] - [Xss; 1] %=0!


%% Using just Bs in augmented state space, not Bs*u:

% %     In = eye(size(As,1) + size(Bs,2), size(As,2) + size(Bs,1));
% %     EA = In;
% %     Atil = zeros(size(EA,1), size(EA,2), n);
% %     for i = n:-1:1
% %         Atil(:,:,i) = [As(:,:,i), Bs(:,:,i); zeros(size(Bs,2), size(EA,1))];
% %         EA = expm(Atil(:,:,i)*ts(i))*EA;
% %     end
% % 
% % 
% %     IHC = eye(size(EA));
% %     IHC(1:ns,1:ns) = obj.IHC;
% % 
% %     depStates = find(sum(abs(eye(ns)-obj.Is(:,:,1)),2));
% %     nd = length(depStates);
% %     EA(depStates, :, :) = [];
% %     EA(:,depStates,:) = [];
% %     IHC(depStates, :, :) = [];
% %     IHC(:,depStates,:) = [];
% % 
% %     In = eye(size(IHC));
% %     Zn = null(In-IHC*EA);
% % 
% %     k = Zn(size(As,1)-nd+1:end,:)\u(:,:,1);
% %     Xss = Zn(1:end -size(Bs,2),:)*k
% % 
% %     EA(1:end-size(Bs,2),1:end-size(Bs,2))*Xss - Xss %=0!




    %% from steady-state solution, go through and find states at each subinterval
    Xs(:,1) = [Xss; 1];
    for i=1:n
        Xs(:,i+1) = expAtil(:,:,i)*Xs(:,i);
    end
    Xs=Xs(1:end-1,:);

    if sum(isnan(Xss))
        ME = MException('resultisNaN:noSuchVariable', ...
                           'Waveform steady-state solution resulted in NaN');
        throw(ME);
    end

    if all(dts == 0)
        % only save the result if it is being evaluated without
        % perturbation
        obj.Xs = Xs;
    end


end

