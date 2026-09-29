function [J, dt] = analyticalJacobian(obj)
%analyticalJacobian calculate discrete Jacobian analytically using exact cyclic steady-state sensitivities
%   
%   [J] = analyticalJacobian(obj)
%   for the SMPSim object obj, calculates discrete Jacobian J.  The
%   Jacobian is a 3-dimensional matrix with structure J(i, j, k): partial
%   derivative of X(i,j) due to perturbations in t(k). 
%
%   See Also SMPSim.discreteJacobian, SMPSim.findValidSteadyState

    assert(~isempty(obj.Xs), 'steady-state solution must be found before computing jacobian');
    
    nt = size(obj.Xs,2)-1;
    ns = size(obj.Xs,1);
    
    J = zeros(ns, nt, nt);
    dt = obj.ts; 
    
    As = obj.As;
    Bs = obj.Bs;
    ts = obj.ts;
    u = obj.fullu;
    
    Phi = zeros(ns, ns, nt);
    X_cont = zeros(ns, nt); % Store continuous states at the end of each interval
    
    %% Precompute the State Transition Matrices (Phi_k) and Continuous States
    if ~all(diff(obj.Is,1,3),'all')
        % method 1 (fast method, consistent Is)
        for k = 1:nt
            At = [As(:,:,k), Bs(:,:,k)*u(:,:,k); zeros(1, ns+1)] * ts(k);
            eatilde = expm(At);
            eat = eatilde(1:end-1, 1:end-1);
            convInt = eatilde(1:end-1, end);
            
            X_cont(:, k) = eat * obj.Xs(:, k) + convInt;
            Phi(:,:,k) = obj.Is(:,:,k) * eat;
        end
    else
        % method 2 (slower method, varying Is)
        for k = 1:nt
            dS = sum(abs(eye(ns)-obj.Is(:,:,k)),2) ~= 0;
            At = [As(~dS,~dS,k), Bs(~dS,:,k)*u(:,:,k); zeros(1, sum(~dS)+1)] * ts(k);
            eatilde = expm(At);
            eat = eatilde(1:end-1, 1:end-1);
            convInt = eatilde(1:end-1, end);
            
            fullPhi = zeros(ns, ns);
            fullPhi(~dS,~dS) = eat;
            fullPhi(dS,:) = obj.Is(dS,:,k);
            Phi(:,:,k) = fullPhi;
            
            X_prev = obj.Xs(:, k);
            X_cont(~dS, k) = eat * X_prev(~dS) + convInt;
            X_cont(dS, k) = 0; 
        end
    end
    
    %% Compute the Global Cycle Inverse
    Phi_total = eye(ns);
    for k = 1:nt
        Phi_total = Phi(:,:,k) * Phi_total;
    end
    
    % Use robust pseudo-inverse
    inv_cyclic = pinv(eye(ns) - Phi_total); 
    
    %% Build the Jacobian
    for i = 1:nt
        % Local perturbation at the end of perturbed interval i
        % Derivative of continuous state before the jump
        x_dot_cont = As(:,:,i) * X_cont(:, i) + Bs(:,:,i) * u(:,:,i);
        
        % Apply the jump mapping to the derivative
        x_dot_i = obj.Is(:,:,i) * x_dot_cont;
        
        % Propagate to the cycle boundary (x_0)
        v = x_dot_i;
        for k = i+1:nt
            v = Phi(:,:,k) * v;
        end
        dx0_dti = inv_cyclic * v;
        
        % Populate the Jacobian for all boundaries
        w = dx0_dti;
        for k = 1:nt
            % Store state derivative at BEGINNING of interval k to match discreteJacobian format
            J(:, k, i) = w; 
            if k == i
                w = Phi(:,:,k) * w + x_dot_i;
            else
                w = Phi(:,:,k) * w;
            end
        end
    end
end
