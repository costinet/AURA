function [ ys, t ] = SS_WF_Reconstruct( Xss, As, Bs, ts, u )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    Ts = sum(ts);
    t = linspace(0,Ts,10000);
    totalts = [0 cumsum(ts)];
    dt = t(2)-t(1);

    nsub = size(As,3);
    n = size(As,1);

    ys = zeros(n, length(t));

    ys = Xss(:,1);


    for i = 1:nsub %for each subinterval
        tmin = totalts(i);
        tmax = totalts(i+1);

        ti = find(t<tmax & t>=tmin);
        if(ti)
            if(length(ti)>1)
                SS = ss(As(:,:,i), Bs(:,:,i), ones(1,n), 0);
%                 [y, tx, x] = lsim(SS, u*ones(size(ti)), t(ti)-t(ti(1)), ys(:,max(ti(1)-1,1)));
                [~, ~, x] = lsim(SS, u*ones(size(ti)), t(ti)-t(ti(1)), ys(:,end));
                ys(:,ti) = x';

            end
        end
    end

    ys = [Xss(:,1) ys];
    
    if sum(isnan(ys))
        ME = MException('resultisNaN:noSuchVariable', ...
                           'Waveform reconstruction resulted in NaN');
        throw(ME);
    end

end
