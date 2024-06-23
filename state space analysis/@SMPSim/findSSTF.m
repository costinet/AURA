function Gz = findSSTF(obj, tp, oi, tm, ut)
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
%   Gz = findSSTF(obj, Km) instead uses a vector Km of gains from each
%   switching interface, so that multiple simultaneous perturbations can be
%   considered.  numel(Km) == numel(obj.ts).
%
%   Gz = findSSTF(___, oi) only outputs the state signals referenced by
%   indices in oi.
%
%   Gz = findSSTF(__, tm) also includes a measurement interval tm
%   at which all state outputs are measured.  If tm is ommitted,
%   measurements are taken after the final interval, tm=length(obj.ts)
%
%   see also SMPSim, ss   

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
%     if size(tp) == size(obj.converter.fullts,2)
%         tps = tp;
%     elseif numel(tp) == 1
%         tps = zeros(1,size(obj.converter.fullts,2));
%         tps(tp) = 1;
%     tps = find(tps);

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
    ActualStopTime = find(ints == tm, 1, 'last');
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

    nt = size(As,3);
    ns = size(As,1);
    ni = size(Bs,2);

    % In = eye(ns + ni, ns + ni);
    % EA = In;
    % Atil = zeros(size(EA,1), size(EA,2), ActualStopTime);
    % for i = ActualStopTime:-1:1
    %     Atil(:,:,i) = [As(:,:,i), Bs(:,:,i); zeros(ni, size(EA,1))];
    %     expAtil(:,:,i) = expm(Atil(:,:,i)*ts(i));
    %     EA = expAtil(:,:,i)*EA;
    % end
    % 
    % PHI = EA(1:ns,1:ns);

    % PHI = eye(ns);
    % for i = ActualStopTime:-1:1
    %     PHI = PHI*expm(As(:,:,i)*ts(i));
    % end
    % 
    % PHI1 = PHI;
    % 
    % PHI = eye(ns);
    % for i = ActualStopTime:-1:1
    %     depStates = sum(abs(eye(ns)-obj.Is(:,:,i)),2)  ~= 0;
    %     EA = expm(As(~depStates,~depStates,i)*ts(i));
    %     EA2 = zeros(ns,ns);
    %     EA2(~depStates,~depStates) = EA;
    %     EA2 = obj.Is(:,:,i)*EA2;
    %     PHI = PHI*EA2;
    % end
    % Handle Dependent States?
    

    
    xphat = zeros(ns,1,nt);
    for i = 1:numel(actualPertTime)
        if actualPertTime(i) ~= 0
            XP = obj.Xs(:,actualPertTime(i)+1);
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


    PHI = eye(ns);
    for i = ActualStopTime:-1:1
        if any(i == actualPertTime)
            % every interval after this one already multiplied in
            Aprop(:,:,i) = PHI;
        end
        EA = zeros(ns,ns);
        depStates = sum(abs(eye(ns)-obj.Is(:,:,i)),2)  ~= 0;
        EA(~depStates,~depStates) = expm(As(~depStates,~depStates,i)*ts(i));
        EA= obj.Is(:,:,i)*EA;

        PHI = PHI*EA;

    end


    % Aprop = eye(ns);
    % for i = ActualStopTime:-1:actualTime+1
    %     Aprop = Aprop*expm(As(:,:,i)*ts(i));
    % end
    % % 
    % Aprop = eye(ns);
    % for i = ActualStopTime:-1:actualTime+1
    %     depStates = sum(abs(eye(ns)-obj.Is(:,:,i)),2)  ~= 0;
    %     EA = expm(As(~depStates,~depStates,i)*ts(i));
    %     EA2 = zeros(ns,ns);
    %     EA2(~depStates,~depStates) = EA;
    %     EA2 = obj.Is(:,:,i)*EA2;
    %     Aprop = Aprop*EA2;
    % end


    if any(ActualStopTime > actualPertTime(actualPertTime>0))
        GAMMA = sum(pagemtimes(Aprop,xphat(:,:,1:size(Aprop,3))),3);
    else
        GAMMA = zeros(ns,1);
    end

    Ts = sum(obj.ts(1:ActualStopTime));

    C = zeros(ns, ns);
    C(sub2ind(size(C),oi,oi)) = 1;
    C(all(C==0,2),:) = [];
    Gz = ss(PHI, GAMMA, C, zeros(size(C,1), 1), Ts);

end