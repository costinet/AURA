function [ xs, t, ys, xd ] = simTransient(obj, nPeriods, Xss, tsteps)
% simTransient step transient waveform reconstruction for periodic switched systems
%
% [ xs, t, ys, xd ] = simTransient(obj, nPeriods) produces
% time domain system states xs, outputs ys, and time vector t for the
% switched system described by the SMPSim object obj.  nPeriods specifies
% the number of periods to carry out the simulation.  obj is a SMPSim
% object, whose internal steady-state timing vector and states will be
% used.
%
% simTransient(obj, nPeriods, Xss) Specifying Xss initializes the states to
% vector Xss, rather than the interally-stores state vector.  
% 
% simTransient(obj, nPeriods, Xss, tsteps) the additional tsteps
% argument allows the user to specify the number of linearl-spaced points 
% in the timing vector, t.

    if nargin < 3 || isempty(Xss) 
        Xss = obj.Xs;
    else
        Xss_stored = obj.Xs;
        obj.Xs = Xss;
    end
    
    if nargin < 4 
        tsteps = min(max(sum(obj.ts)/min(obj.ts)*100, 10e3),100e3);
        tsteps = round(tsteps);
    end
    
    xs = zeros(nPeriods*tsteps,size(Xss,1));
    ys = zeros(nPeriods*tsteps,size(obj.Cs,1));
    t = zeros(nPeriods*tsteps,1);

    xd = zeros(size(Xss,1),nPeriods);

    n = size(obj.As,3);
    ns = size(obj.Bs,1);
    In = eye(ns + 1);
    u=obj.fullu;
    
    % Calculate discrete State Transition Matrix (Augmented)
    PHIaug = In;
    for i = n:-1:1
       Atil = [obj.As(:,:,i), obj.Bs(:,:,i)*u(:,:,i); zeros(1, size(PHIaug,1))];
        expAtil = expm(Atil*obj.ts(i));
        PHIaug = PHIaug*expAtil;
    end

    
    endPrd = 1;
    for i = 1:nPeriods
        [ xsT, tT, ysT ] = simSinglePeriod(obj, Xss, tsteps);
        xs(endPrd:endPrd+size(xsT,2)-1,:) = xsT';
        t(endPrd:endPrd+numel(tT)-1,:) = tT' + t(max(1,endPrd-1));
        ys(endPrd:endPrd+size(ysT,2)-1,:) = ysT';
    
        endPrd = endPrd+numel(tT);
        xd(:,i) = Xss;
        % Xss = xs(endPrd-1,:);

        % Calculate discrete per-period sample points explicitly to
        % minimize cumulative error across subinterval transitions/periods.
        Xss = PHIaug*[Xss; 1];
        Xss = Xss(1:end-1);
    end

    if exist('Xss_stored', 'var')
        obj.Xs = Xss_stored;
    end
end

function [ xs, t, ys ] = simSinglePeriod(obj, X0, tsteps)
% waveform reconstruction for periodic switched systems



    
    As = obj.As;
    Bs = obj.Bs;
    ts = obj.ts;
    % u = obj.u;
    Cs = obj.Cs;
    Ds = obj.Ds;
    u = obj.fullu;
    
    try
        %wrapped in try block for debugging purposes
        Ts = sum(ts);
        totalts(1) = 0;
        totalts(2:length(ts)+1) = cumsum(ts);
    
        nsub = size(As,3);
        n = size(As,1);
        
        if nargin == 1
            tsteps = min(max(Ts/min(ts)*100, 10e3),100e3);
    
        end
    
        t = linspace(0,Ts,tsteps);
    
        xs = zeros(n, length(t));
        xs(:,1) = X0;
        
        ys = zeros(size(Cs,1), length(t));
    
        for i = 1:nsub %for each subinterval
            tmin = totalts(i);
            tmax = totalts(i+1);
    
            if i < nsub
                ti = find(t<tmax & t>=tmin);
            elseif i == nsub
                ti = find(t<=tmax & t>=tmin);
            end
    
            if length(ti) > 1
                SS = ss(As(:,:,i), Bs(:,:,i), Cs(:,:,i), Ds(:,:,i));
                [y, ~, x] = lsim(SS, u(:,:,i)*ones(size(ti)), t(ti)-t(ti(1)), X0);
    
                xs(:,ti) = x';
                ys(:,ti) = y';
    
                X0 = x(end,:)';
    
            elseif isscalar(ti)% special case if only one timestep in range
               tinter = linspace(0, t(2), 100);
               SS = ss(As(:,:,i), Bs(:,:,i), Cs(:,:,i), Ds(:,:,i));
               [y, ~, x] = lsim(SS, u(:,:,i)*ones(size(tinter)), tinter, xs(:,max(ti(1)-1,1)));
               ind = find(abs(tinter + tmin - t(ti)) == min(abs(tinter + tmin - t(ti))),1);
               xs(:,ti) =x(ind,:)';
               ys(:,ti) = y(ind,:)';
    
            end
        end
    
        xs = [X0 xs(:,1:end-1)];
        
        if sum(isnan(xs))
            ME = MException('resultisNaN:noSuchVariable', ...
                               'Waveform reconstruction resulted in NaN');
            throw(ME);
        end
    catch ERROR
        rethrow(ERROR)
    end

end
