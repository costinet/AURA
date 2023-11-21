function [ xs, t, ys ] = SS_WF_Reconstruct(obj, tsteps)
% Steady-state waveform reconstruction for periodic switched systems
%
% [ xs, t, ys ] = SS_WF_Reconstruct(obj) produces
% time domain system states xs, outputs ys, and time vector t for the
% switched system described by the SMPSim object obj
%
% [ xs, t, ys ] = SS_WF_Reconstruct(obj, tsteps) the additional tsteps
% argument allows the user to specify the number of linearl-spaced points 
% in the timing vector, t.


Xss = obj.Xs;
if isempty(Xss)
    obj.steadyState;
    Xss = obj.Xs;
end
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
%         [~,dt] = obj.converter.getDeltaT;
%         tsteps = ts./dt
    end
%     if nargin < 7
%         Ds = zeros(1,1,nsub);
%     end
%     if nargin < 6
%         Cs = ones(1,n, nsub);
%     end
    
    t = linspace(0,Ts,tsteps);

    xs = zeros(n, length(t));
    xs(:,1) = Xss(:,1);
    
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
            [y, ~, x] = lsim(SS, u(:,:,i)*ones(size(ti)), t(ti)-t(ti(1)), Xss(:,i));

            xs(:,ti) = x';
            ys(:,ti) = y';

        elseif length(ti) == 1 % special case if only one timestep in range
           tinter = linspace(0, t(2), 100);
           SS = ss(As(:,:,i), Bs(:,:,i), Cs(:,:,i), Ds(:,:,i));
           [y, ~, x] = lsim(SS, u(:,:,i)*ones(size(tinter)), tinter, xs(:,max(ti(1)-1,1)));
           ind = find(abs(tinter + tmin - t(ti)) == min(abs(tinter + tmin - t(ti))),1);
           xs(:,ti) =x(ind,:)';
           ys(:,ti) = y(ind,:)';

        end
    end

    xs = [Xss(:,1) xs(:,1:end-1)];
    
    if sum(isnan(xs))
        ME = MException('resultisNaN:noSuchVariable', ...
                           'Waveform reconstruction resulted in NaN');
        throw(ME);
    end
catch ERROR
    rethrow(ERROR)
end

end
