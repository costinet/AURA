function [ Xs] = SS_Soln(obj, Xi, Bi)
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


    As = obj.As;
    Bs = obj.Bs;
    ts = obj.ts;
    u = obj.u;

    if(nargin == 1)
        tryOpt = 0;
    else
        tryOpt = obj.tryOpt;
    end

    n = size(As,3);
    ns = size(Bs,1);

    condAs= zeros(1,n);

    % Compute expAS = e^(A_i*t_i)
    % where A_i is A matrix for state i period
    % t_i is the length of time in state i in period
    % states are in order i = 1,2,3...
    expAs = zeros(ns,ns,n);
    for i=1:n
        condAs(i) = cond(As(:,:,i)); % check the condition of inversion for A matrix
        expAs(:,:,i) = expm(As(:,:,i)*ts(i)); % e^(A*t)
    end

    
    % Compute cumProdExp = e^(A_(i) * t_(i)) * e^(A_(i-1) * t_(i-1))
    % ... e^(A_(1) * t_(1)) * I
    
    % I is the identity matrix
    cumProdExp = zeros(ns,ns,n);
    for i=1:n
        cumProdExp(:,:,i) = eye(ns);
        for k = 1:i
            cumProdExp(:,:,i) = expAs(:,:,k)*cumProdExp(:,:,i);
        end
    end

    
    % Compute cumProdExpRev = e^(A_(n) * t_(n)) * e^(A_(i+2) * t_(i+2))
    % ... e^(A_(i+1) * t_(i+1)) * I
    % Unless i = n then cumProdExpRev = I

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

    fresp = zeros(ns,1,n); % Forced response: A^-1(expm(A*T)-eye)*B*u
        for i=1:n
           [frespNew, intEAt] = obj.forcedResponse(As(:,:,i), expAs(:,:,i), Bs(:,:,i), u, ts(i), 1);
           fresp(:,:,i) = frespNew;
           obj.oldAs(:,:,i) = As(:,:,i);
           obj.oldIntEAt(:,:,i) = intEAt;
           obj.oldts(i) = ts(i);
           RHSsum = RHSsum +  cumProdExpRev(:,:,i)*frespNew;
        end    
    hi = -6; % Intial guess (might be able to educated guess this)
    
    % This effectivly places a large shut resistor on all caps and series
    % resistors on all all inductors:
    if cond(eye(ns) - cumProdExp(:,:,n))>1*10^9 % 10^-9 is an educated guess 
        Xss = ((1-10^(hi))*eye(ns) - cumProdExp(:,:,n))^-1*RHSsum; % Calculate the educated guess value

        %%%% Use cond() to see where how 'able' the matrix is to converge %%%%
        %%% Can also try and use optimization commented out below %%%%
        
        % Loop through increaseing hi until the condition of the inverted
        % matrix reaches a 'large' value
        while cond((1-10^(hi))*eye(ns) - cumProdExp(:,:,n))<1*10^9
            hi = hi-1;
            Xss = ((1-10^(hi))*eye(ns) - cumProdExp(:,:,n))^-1*RHSsum;
        end
    else
        Xss = (eye(ns) - cumProdExp(:,:,n))^-1*RHSsum;
    end
        
    %% Hacky solution when result is off because (eye(ns) - cumProdExp(:,:,n)) is non-invertable
    % Use optimization to solve without inversion through error minimization
    % --> Still needs tweaking
    %{
if(tryopt)
        if(abs(Xss(1) - Xi(1)) > 2)
            warning(['Unable to find SS solution directly for Xi = ' num2str(Xi)]);
            Xss = Xi;
            options = optimoptions('lsqlin','algorithm','trust-region-reflective','Display','none');
            A = (eye(ns) - cumProdExp(:,:,n));
            b = RHSsum;
            Xss = lsqlin(A, b, [],[],[],[], Xi.*Bi, Xi.*(2-Bi), Xss, options);
        end
    end
%}
    %% from steady-state solution, go through and find states at each subinterval
    Xs(:,1) = Xss;
    for i=1:n
        Xs(:,i+1) = expAs(:,:,i)*Xs(:,i) + fresp(:,:,i);
    end

    if sum(isnan(Xss))
        ME = MException('resultisNaN:noSuchVariable', ...
                           'Waveform steady-state solution resulted in NaN');
        throw(ME);
    end

    obj.Xs = Xs;


end

