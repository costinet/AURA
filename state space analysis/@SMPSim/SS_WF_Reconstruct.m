function [ xs, t, ys ] = SS_WF_Reconstruct(obj, tsteps)
% Steady-state waveform reconstruction for periodic switched systems
%
% [ xs, t, ys ] = SS_WF_Reconstruct( Xss, As, Bs, ts, u, Cs, Ds ) produces
% time domain system states xs, outputs ys, and time vector t for the
% switched system described by 
%
% dx/dt = Ai*x(t) + Bi*u ,
% y(t) = Ci*x(t) + Di*u ,
%
% for the ith interval.  The function repetatively applies lsim to simulate
% the linear system within each subinterval, while matching boundary
% conditions between each.
%
% Xss is a 2D matrix with ns rows and n+1 columns, where ns is the
% number of states in the system and n is the number of switching
% intervals.  The first and last column of Xs should be identical,
% corresponding to a valid steady-state solution.
% As is a 3-dimensional matrix of values for Ai, where As(:,:,i) is the 2D
% square matrix Ai during the ith interval.  
% Bs is a 3-dimensional matrix of values for Bi, where Bs(:,:,i) is the 2D
% matrix/vector Bi during the ith interval.
% Cs is a 3-dimensional matrix of values for Ci, where Cs(:,:,i) is the 2D
% square matrix Ci during the ith interval.  
% Ds is a 3-dimensional matrix of values for Di, where Ds(:,:,i) is the 2D
% vector Di during the ith interval.
% ts is a vector of the time durations of each inverval
% u is the (assumed constant) independent input vector
% tsteps is the number of timesteps in one period, which is auto-selected
% if left blank, to be between 10k and 100k while trying to keep at least
% 100 steps per switching subinterval

Xss = obj.Xs;
As = obj.As;
Bs = obj.Bs;
ts = obj.ts;
u = obj.u;
Cs = obj.Cs;
Ds = obj.Ds;

try
    %wrapped in try block for debugging purposes
    Ts = sum(ts);
    totalts = [0 cumsum(ts)];

    nsub = size(As,3);
    n = size(As,1);
    
    if nargin == 1
        tsteps = min(max(Ts/min(ts)*100, 10e3),100e3);
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
            [y, ~, x] = lsim(SS, u*ones(size(ti)), t(ti)-t(ti(1)), Xss(:,i));

            xs(:,ti) = x';
            ys(:,ti) = y';

        elseif length(ti) == 1 % special case if only one timestep in range
           tinter = linspace(0, t(2), 100);
           SS = ss(As(:,:,i), Bs(:,:,i), Cs(:,:,i), Ds(:,:,i));
           [y, ~, x] = lsim(SS, u*ones(size(tinter)), tinter, xs(:,max(ti(1)-1,1)));
           ind = find(abs(tinter + tmin - t(ti)) == min(abs(tinter + tmin - t(ti))),1);
           xs(:,ti) = x(ind);
           ys(:,ti) = y(ind);

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
