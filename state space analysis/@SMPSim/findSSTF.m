function Gz = findSSTF(obj, tp, tm, ut)
%findSSTF finds the small-signal transfer function from time interval tp to
%all states
%   Gz = findSSTF(obj, tp) finds the z-domain transfer function from a time
%   interval specified by tp to all states of the system. obj is an object 
%   of type SMPSim for which a steady-state solution has been previously
%   found.  tp is an integer index into the originally-specified switching
%   times of the converter.
%   
%   Gz is a discrete time state space model containing ns z-domain transfer
%   functions, where the ith transfer function is the small signal gain
%   from time interval tp to state i.
%
%   Gz = findSSTF(obj, tp, tm) also includes a measurement interval tm
%   at which all state outputs are measured.  If tm is ommiteed,
%   measurements are taken after the final interval, tm=length(obj.ts)
%
%   see also SMPSim, ss   

    if nargin==2
        tm = size(obj.converter.fullts,2);
    else
         assert( numel(tm) == 1 && tm > tp && tm<=size(obj.converter.fullts,2), 'tm must be an integer index into the controlled switch times after tp')
    end

    if isempty(tp)
        Gz = 0;
        return
    end
    assert( numel(tp) == 1 && tp > 0 && tp<size(obj.converter.fullts,2), 'tp must be an integer index into the controlled switch times')
%     if size(tp) == size(obj.converter.fullts,2)
%         tps = tp;
%     elseif numel(tp) == 1
%         tps = zeros(1,size(obj.converter.fullts,2));
%         tps(tp) = 1;
%     tps = find(tps);

    if isempty(obj.Xs)
        obj.steadyState;
    end

    if nargin < 4
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
    ActualStopTime = find(ints == tm, 1, 'last');

    nt = size(As,3);
    ns = size(As,1);
    ni = size(Bs,2);

    In = eye(ns + ni, ns + ni);
    EA = In;
    Atil = zeros(size(EA,1), size(EA,2), ActualStopTime);
    for i = ActualStopTime:-1:1
        Atil(:,:,i) = [As(:,:,i), Bs(:,:,i); zeros(ni, size(EA,1))];
        expAtil(:,:,i) = expm(Atil(:,:,i)*ts(i));
        EA = EA*expAtil(:,:,i);
    end

    PHI = EA(1:ns,1:ns);

    
    actualTime = find(ints == tp, 1, 'last');

    XP = obj.Xs(:,actualTime+1);
    xphat = As(:,:,actualTime)*XP + Bs(:,:,actualTime)*u(:,:,actualTime) - ...
        As(:,:,actualTime+1)*XP + Bs(:,:,actualTime+1)*u(:,:,actualTime+1);

    Aprop = eye(ns);
    for i = ActualStopTime:-1:actualTime
        Aprop = expm(As(:,:,i)*ts(i))*Aprop;
    end

    GAMMA = Aprop*xphat;

    Ts = sum(obj.ts(1:ActualStopTime));

    Gz = ss(PHI, GAMMA, eye(ns), zeros(ns, 1), Ts);

end