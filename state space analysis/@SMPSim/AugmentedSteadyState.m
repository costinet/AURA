function [ Xs] = AugmentedSteadyState(obj, dts)
%Steady-state solution of switched system using state-space matrices
%
%   [ Xs] = AugmentedSteadyState(obj) finds the state values Xs in
%   steady-state for the switched system described by the SMPSim object obj
%
%   [ Xs] = AugmentedSteadyState(obj,dts) finds the steady-state solution
%   with perturbations dts to the timing of each interval.
%
%   The result, Xs is a 2D matrix with ns rows and n+1 columns, where ns is the
%   number of states in the system and n is the number of switching
%   intervals.  The first and last column of Xs should be identical,
%   corresponding to a valid steady-state solution.

   
    assert(isa(obj,'SMPSim'),'AugmentedSteadyState requires an object of type SMPSim as its first argument');

    if nargin <2
        dts = obj.ts*0;
    else
        assert(length(dts) == length(obj.ts), 'perturbations dts must be the same length as obj.ts');
    end

    As = obj.As;
    Bs = obj.Bs;
    ts = obj.ts + dts;

    if isempty(obj.IHC)
        obj.IHC = eye(size(As,1));
    end
    
    u = obj.fullu;

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


%     PHI = In(~depStates,~depStates)-IHC(~depStates,~depStates)*EA(~depStates,~depStates);
%     Zn = null(PHI(1:end-1,1:end-1));
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

