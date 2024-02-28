function [ avgXs, avgYs, Ints ] = ssAvgs(obj, Xss)
%ssAvgs Compute average values of states and outputs in steady-state
%   
%   [ avgXs, avgYs, Ints ] = ssAvgs(obj)
%   for SMPSim object, computes the vector of average states avgXs and
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
    %           ns      1       ns      nc
    %           __      __      __      __  
    %   A = [   As      Bs*u    0       0       |ns
    %           0       0       0       0       |1    
    %           I       0       0       0       |ns
    %           Cs      Ds*u    0       0]      |nc

    Atil = [As, pagemtimes(Bs,u), zeros(ns,ns+nc,n);
            zeros(1,2*ns+1+nc,n);
            repmat(eye(ns),[1 1 n]), zeros(ns,ns+1+nc,n);
            Cs, pagemtimes(Ds,u), zeros(nc,nc+ns,n)];
    Xtil = [Xss; ones(1,n+1); zeros(ns,n+1); zeros(nc,n+1)];

    Ints = zeros(ns*2+nc+1,n);
    for i = 1:n
        Ints(:,i) =  expm(Atil(:,:,i)*ts(i))*Xtil(:,i);
    end

    avgXs = sum(Ints(ns+2:2*ns+1,:),2)/sum(ts);
    avgYs =  sum(Ints(2*ns+2:end,:),2)/sum(ts);   
    Ints = Ints(ns+2:end,:);

end

