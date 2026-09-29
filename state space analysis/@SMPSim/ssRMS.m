function [rmsXs, rmsYs] = ssRMS(obj, Xss)
%ssRMS Compute rms values of states and outputs in steady-state
%   
%   [ rmsXs, rmsYs, Ints ] = ssRMS(obj)
%   for SMPSim object, computes the vector of rms states rmsXs and
%   outputs avgYs using the current, stored steady-state solution.
%   Ints gives the integrals of each state and output over each time
%   interval in a (Ns+Ni x Nt) matrix for a system with Ns states, Ni
%   inputs, and Nt time intervals
%   
%   [ ... ] = ssAvgs(obj, Xss)
%   uses instead the supplied steady-state initial vector Xss
%
%   See Also SMPSim.steadyState, SMPSim.ssAvgsInvert

    As = obj.As;
    Bs = obj.Bs;
    Cs = obj.Cs;
    Ds = obj.Ds;
    Is = obj.Is;
    ts = obj.ts;
    u = obj.fullu;

    n = size(As,3);
    ns = size(Bs,1);
    nc = size(Cs,1);

    if nargin<2
        Xss = obj.Xs;
        if isempty(obj.Xs) || size(Xss,2) ~= n 
            Xss = obj.steadyState;
        end
    end


    
    %% Augmented State matrix
    %           ns      1       ns^2    ns      nc
    %           __      __      __      __      __
    %   A = [   As      Bs*u    0       0       0   |ns
    %           0       0       0       0       0   |1    
    %           Bq              Aq                  |ns^2
    %                           Sx                  |ns
    %           Syx     Syd     Syq                 |nc
                
    %        ns 1  ns^2 ns    nc
    %   x = [x, 1, q,   Ix^2, Iy^2]

    Sx = zeros(ns,ns^2);
    for k =1:n+1
        if k<=n
            Aq(:,:,k) = kron(eye(ns),As(:,:,k)) + kron(As(:,:,k),eye(ns));
            Bq(:,:,k) = kron(eye(ns),Bs(:,:,k)*u(:,:,k)) + kron(Bs(:,:,k)*u(:,:,k),eye(ns));
    
    
            Syq = zeros(nc,ns^2);
            % for i = 1:nc
            %     ci = Cs(i,:,k);
            %     Syq(i,:,k) = kron(ci,ci);
            %     % ci = Cs(i,:,k) * Is(:,:,k);
            %     % Syq(i,:,k) = kron(ci,ci);
            % end
            for i = 1:nc
                ci = Cs(i,:,k) * Is(:,:,k);

                Qi = ci.' * ci;
                Syq(i,:,k) = Qi(:).';
            end
    
        
            % Syx(:,:,k) = 2*diag(Ds(:,:,k)*u(:,:,k))*Cs(:,:,k);
            Syx(:,:,k) = 2 * diag(Ds(:,:,k)*u(:,:,k)) ...
                * Cs(:,:,k) * Is(:,:,k);
            Syd(:,:,k) = Ds(:,:,k)*u(:,:,k).^2;

            for i=1:ns
                % Sx(i,:,k) = kron(Is(i,:,k),Is(i,:,k));
                ri = Is(i,:,k);
                Sx(i,:,k) = reshape(ri.' * ri, 1, []);
            end
        end

        qt = Xss(:,k)*Xss(:,k)';
        q(:,k) = qt(:);


    end

    % Sx = zeros(ns,ns^2);
    % for i=1:ns
    %     % Sx(i,(i-1)*ns+i) = 1;
    %     Sx(i,:,k) = kron(Is(i,:,k),Is(i,:,k));
    % end
    % Sx = repmat(Sx,[1 1 n]);


    % Original, rms-free augmentation
    % Atil = [As, pagemtimes(Bs,u), zeros(ns,ns+nc,n);
    %         zeros(1,2*ns+1+nc,n);
    %         repmat(eye(ns),[1 1 n]), zeros(ns,ns+1+nc,n);
    %         Cs, pagemtimes(Ds,u), zeros(nc,nc+ns,n)];
    % Xtil = [Xss; ones(1,n+1); zeros(ns,n+1); zeros(nc,n+1)];
       

    %
    % Ap = [Atil, zeros(ns*2+1+nc, ns^2+ns+nc, n);
    %     Bq, zeros(ns^2,1+ns+nc,n), Aq, zeros(ns^2,ns+nc,n);
    %     zeros(ns,2*ns+1+nc,n), Sx, zeros(ns, ns+nc,n);
    %     Syx, Syd, zeros(nc,ns+nc,n), Syq, zeros(nc, ns+nc,n)];
    % Xp = [Xss; ones(1,n+1); zeros(ns,n+1); zeros(nc,n+1); q; zeros(ns+nc,n+1)];

    Abase = [As, pagemtimes(Bs,u);
        zeros(1,ns+1,n)];

    Ap = [Abase, zeros(ns+1, ns^2+ns+nc, n);
        Bq, zeros(ns^2,1,n), Aq, zeros(ns^2,ns+nc,n);
        zeros(ns,ns+1,n), Sx, zeros(ns, ns+nc,n);
        Syx, Syd, Syq, zeros(nc, ns+nc,n)];
    Xp = [Xss; ones(1,n+1); q; zeros(ns+nc,n+1)];

    Ints = zeros(ns*2+nc+ns^2+1,n);
    for i = 1:n
        Ints(:,i) =  expm(Ap(:,:,i)*ts(i))*Xp(:,i);
    end

    % avgXs = sum(Ints(ns+2:2*ns+1,:),2)/sum(ts);
    % avgYs =  sum(Ints(2*ns+2:2*ns+1+nc,:),2)/sum(ts);   
    rmsXs = sqrt(sum(Ints(ns+ns^2+2:2*ns+ns^2+1,:),2)/sum(ts));
    rmsYs = sqrt(sum(Ints(2*ns+ns^2+2:2*ns+ns^2+1+nc,:),2)/sum(ts));
    % Ints = Ints(ns+2:2*ns+1+nc,:);

end

