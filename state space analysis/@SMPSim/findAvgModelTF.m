function Gs = findAvgModelTF(obj, tp, oi)
%findAvgModelTF finds the small-signal transfer function from time interval tp to
%all states using state space averaging.  
% For discrete time transfer functions, use findSSTF()
%
%   Gz = findAvgModelTF(obj, tp) finds the z-domain transfer function from a time
%   interval specified by tp to all states of the system. obj is an object 
%   of type SMPSim for which a steady-state solution has been previously
%   found.  tp is an integer index into the originally-specified switching
%   times of the converter.
%   
%   Gz is a average state space model containing ns s-domain transfer
%   functions, where the ith transfer function is the small signal gain
%   from time interval tp to state i.
%
%   Gz = findAvgModelTF(obj, Km) instead uses a vector Km of gains from each
%   switching interface, so that multiple simultaneous perturbations can be
%   considered.  numel(Km) == numel(obj.ts).
%
%   Gz = findAvgModelTF(___, oi) only outputs the state signals referenced by
%   indices in oi.
%
%   see also SMPSim, findSSTF, ss 

    if ~exist('tm', 'var')
        tm = size(obj.converter.fullts,2);
    else
         assert( isscalar(tm) && all(tm > tp) && tm<=size(obj.converter.fullts,2), 'tm must be an integer index into the controlled switch times after tp')
    end

    if isempty(tp)
        Gz = 0;
        return
    end

    if isscalar(tp)
        assert( tp > 0 && tp<size(obj.converter.fullts,2), 'tp must be an integer index into the controlled switch times')
        Km = ones(1,numel(obj.ts));
    else
        Km = tp;
        assert(numel(tp) == size(obj.converter.fullts,2), 'Km must be specified as a vector of gains of the same lenght as the controlled switching times')
    end

    if ~exist('oi', 'var') || isempty(oi)
        oi = 1:size(obj.As,1);
    end

    if isempty(obj.Xs)
        obj.steadyState;
    end

    if nargin <= 4
        ut = 1;
    else 
        error('functionality for input transfer function not implemented')
    end

    As = obj.As;
    Bs = obj.Bs;
    ts = obj.ts;

    if isempty(obj.IHC)
        obj.IHC = eye(size(As,1));
    end
    
    u = obj.fullu;

    [~, ints, ~] = obj.converter.getIntervalts();
     if isscalar(tp)
        actualPertTime = find(ints == tp, 1, 'last');
    else
        actualPertTime = 0*obj.ts;
        for i = 1:numel(tp)
            if tp(i)~=0
                actualPertTime(i) = find(ints == i, 1, 'last');
            end
        end
        
     end

    Ts = sum(obj.ts(1:end));
    A = sum(obj.As.*permute(obj.ts,[1 3 2]),3)/Ts;
    B = sum(obj.Bs.*permute(obj.ts,[1 3 2]),3)/Ts;
    XP = -A\B*obj.u;

    nt = size(As,3);
    ns = size(As,1);
    ni = size(Bs,2);
    
    xphat = zeros(ns,1,nt);
    for i = 1:numel(actualPertTime)
        if actualPertTime(i) ~= 0
            % XP = obj.Xs(:,actualPertTime(i)+1);
            if actualPertTime(i) == size(As,3)
                xphat(:,:,actualPertTime(i)) = As(:,:,actualPertTime(i))*XP + Bs(:,:,actualPertTime(i))*u(:,:,actualPertTime(i)) - ...
                    (As(:,:,1)*XP + Bs(:,:,1)*u(:,:,1));
            else
                xphat(:,:,actualPertTime(i)) = As(:,:,actualPertTime(i))*XP + Bs(:,:,actualPertTime(i))*u(:,:,actualPertTime(i)) - ...
                    (As(:,:,actualPertTime(i)+1)*XP + Bs(:,:,actualPertTime(i)+1)*u(:,:,actualPertTime(i)+1));
            end
            xphat(:,:,actualPertTime(i)) = xphat(:,:,actualPertTime(i))*Km(actualPertTime(i));
        end
    end

 

    
    F = sum(xphat(:,:,1:size(obj.As,3)),3);




    C = zeros(ns, ns);
    C(sub2ind(size(C),oi,oi)) = 1;
    C(all(C==0,2),:) = [];
    Gs = ss(A, F, C, zeros(size(C,1), 1));

end