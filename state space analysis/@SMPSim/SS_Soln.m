function [ Xs] = SS_Soln(obj)
%% Steady-state solution of switched system using state-space matrices
%
%   [ Xs] = SS_Soln(obj) finds the state values Xs in steady-state for the 
%   switched system described by the SMPSim object obj
%
%   Compared to AugmentedSteadyState, SS_Soln tries to solve the converter
%   without using augmented state matrices first, and only augments if
%   A(:,:,i) is ill-conditioned.  This is all taken care of in the calls to
%   obj.forcedResponse
%
%   SS_Soln is not in active use, as it is generally slower and less robust
%   than AugmentedSteadyState
%
%   see also SMPSim, SMPSim.AugmentedSteadyState, SMPSim.forcedResponse




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


    As = obj.As;
    Bs = obj.Bs;
    ts = obj.ts;
    
%     if size(obj.u,3) == 1
%         %% HEY THIS SHOULD BE MOVED ELSEWHERE
%         u = repmat(obj.u,1,1,(length(ts)));
%     elseif size(obj.u,3) == length(ts)
%         u = obj.u;
%     else
%         error('invalid input vector u specificed');
%     end
    u = obj.fullu;

%     if(nargin == 1 || obj.tryOpt == 0)
%         tryOpt = 0;
%     else
%         tryOpt = obj.tryOpt;
%     end

    n = size(As,3);
    ns = size(Bs,1);

    condAs= zeros(1,n);

    expAs = zeros(ns,ns,n);
    for i=1:n
        condAs(i) = cond(As(:,:,i));
        expAs(:,:,i) = expm(As(:,:,i)*ts(i));
    end

    cumProdExp = zeros(ns,ns,n);
    for i=1:n
        cumProdExp(:,:,i) = eye(ns);
        for k = 1:i
            cumProdExp(:,:,i) = expAs(:,:,k)*cumProdExp(:,:,i);
        end
    end

    cumProdExpRev = zeros(ns,ns,n);
    for i=1:n
        cumProdExpRev(:,:,i) = eye(ns);
        if(i<n)
            for k = i+1:n
                cumProdExpRev(:,:,i) = expAs(:,:,k)*cumProdExpRev(:,:,i);
            end
        end
    end

    RHSsum = zeros(ns,1);

    fresp = zeros(ns,1,n); %% Forced response: A^-1(expm(A*T)-eye)*B*u
        for i=1:n
           [frespNew, intEAt] = obj.forcedResponse(As(:,:,i), expAs(:,:,i), Bs(:,:,i), u(:,:,i), ts(i), 0);
           fresp(:,:,i) = frespNew;
           obj.oldAs(:,:,i) = As(:,:,i);
           obj.oldIntEAt(:,:,i) = intEAt;
           obj.oldts(i) = ts(i);
           RHSsum = RHSsum +  cumProdExpRev(:,:,i)*frespNew;
        end

    warning('OFF', 'MATLAB:singularMatrix')
    Xss = (eye(ns) - cumProdExp(:,:,n))^-1*RHSsum;
    
    
    %% Check for depdendent states
    if(sum(isnan(Xss)) || sum(isinf(Xss)))
        A = (eye(ns) - cumProdExp(:,:,n));
        zeroCols = sum(A==0,1) == length(A);
        if sum(zeroCols) > 0
            A2 = A(~zeroCols, ~zeroCols);
            b2 = RHSsum(~zeroCols);
            Xss2 = A2\b2;
            Xss = zeros(length(A),1);
            Xss(~zeroCols) = Xss2;
            if ~(sum(isnan(Xss)) || sum(isinf(Xss)))
                % Solution worked, replace the dependent states that were
                % removed
                Xss = obj.Is(:,:,1)*Xss;
%                 tryOpt = 0;
            end
            
%             %% Debugging:
%             if rcond(A2) < 1e-15
%                 test=1
%             end
        end


        %% Hacky solution when result is off because (eye(ns) - cumProdExp(:,:,n)) is non-invertable
        % Use linprog to solve without inversion through error minimization
        % I'm not sure if this addresses any real issue, i.e. if there is a
        % case where this would work but the above would not.
        
        %% 09/28/21 -- Removed this
        %currently left as commented out including Xi and tryopt in case it
        %is needed anytime later.  Should move to another function, if so.
%         if(tryOpt)
%             if(abs(Xss(1) - Xi(1)) > 2 || isnan(Xss(1)))
%                 warning(['Unable to find SS solution directly for Xi = ' num2str(Xi)]);        
%                 X = linprog( Xss2*0+1, [], [], A2, b2);
%                 Xss = obj.Is(:,:,1)*X;
%             end
%         end
    end
    
    warning('ON', 'MATLAB:singularMatrix')

    %% from steady-state solution, go through and find states at each subinterval
    Xs(:,1) = Xss;
    for i=1:n
        Xs(:,i+1) = expAs(:,:,i)*Xs(:,i) + fresp(:,:,i);
%         Xs(:,i+1) = obj.Is(:,:,i)*Xs(:,i+1);
    end

    if sum(isnan(Xss))
        ME = MException('resultisNaN:noSuchVariable', ...
                           'Waveform steady-state solution resulted in NaN');
        throw(ME);
    end

    obj.Xs = Xs;


end

