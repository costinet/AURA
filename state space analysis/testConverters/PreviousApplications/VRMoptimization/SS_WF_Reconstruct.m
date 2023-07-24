function [ xs, t, ys ] = SS_WF_Reconstruct( Xss, As, Bs, ts, u, Cs, Ds )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
try
    Ts = sum(ts);
    tsteps = min(max(Ts/min(ts)*100, 10e3),100e3);
    t = linspace(0,Ts,tsteps);
    totalts = [0 cumsum(ts)];
    dt = t(2)-t(1);

    nsub = size(As,3);
    n = size(As,1);
    
    if nargin < 7
        Ds = zeros(1,1,nsub);
    end
    if nargin < 6
        Cs = ones(1,n, nsub);
    end

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
            [y, tx, x] = lsim(SS, u*ones(size(ti)), t(ti)-t(ti(1)), Xss(:,i));
            
            xs(:,ti) = x';
            ys(:,ti) = y';

        elseif length(ti) == 1 % special case if only one timestep in range
           tinter = linspace(0, t(2), 100);
           SS = ss(As(:,:,i), Bs(:,:,i), Cs(:,:,i), Ds(:,:,i));
           [y, tx, x] = lsim(SS, u*ones(size(tinter)), tinter, xs(:,max(ti(1)-1,1)));
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
