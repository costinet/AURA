function [ ys, t ] = SS_WF_Reconstruct( Xss, As, Bs, ts, u )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
try
    Ts = sum(ts);
    t = linspace(0,Ts,10000);
    totalts = [0 cumsum(ts)];
    % t = [linspace(totalts(1),totalts(2),10) linspace(totalts(2),totalts(3),50) linspace(totalts(3),totalts(4),10) linspace(totalts(4),totalts(5),50)];
    dt = t(2)-t(1);

    nsub = size(As,3);
    n = size(As,1);

    ys = zeros(n, length(t));

    ys = Xss(:,1);


    for i = 1:nsub %for each subinterval
        tmin = totalts(i);
        tmax = totalts(i+1);
    %     if i ==1
    %         tmin = 0;
    %     else
    %         tmin = sum(ts(1:i-1));
    %     end
        ti = find(t<tmax & t>=tmin);
    %     if(i==1)
    %         ys(1,ti) = Xss(1,2);
    %         ys(2,ti) = (Xss(2,i+1)-Xss(2,i))/(tmax-tmin).*(t(ti)-tmin) + Xss(2,i);
    %         ys(3,ti) = Xss(3,1) + cumtrapz((t(ti)-tmin),ys(2,ti));
    %     elseif(i==3)
    %         ys(1,ti) = 0;
    %         ys(2,ti) = (Xss(2,i+1)-Xss(2,i))/(ts(i+1)-ts(i)).*t(ti) + Xss(2,i);
    %     elseif(ti)
        if(ti)
    %           ys(:,ti) = expm(As(:,:,i)*dt)*ys(:,) + As(:,:,i)^-1*(expm(As(:,:,i)*dt) - eye(n))*Bs(:,:,i)*u;
    %         ys(:,ti) = expm(As(:,:,i).*(t(ti)-tmin))*ys(find(t>totalts(i),1) +  As(:,:,i)^-1*(expm(As(:,:,i).*(t(ti)-tmin))) - eye(n))*Bs(:,:,i)*u;

            SS = ss(As(:,:,i), Bs(:,:,i), ones(1,n), 0);
            [y, tx, x] = lsim(SS, u*ones(size(ti)), t(ti)-t(ti(1)), ys(:,max(ti(1)-1,1)));
            ys(:,ti) = x';
    % % %         for j = ti
    % % % % %             if (j<length(t))
    % % % % %                 dt = t(j+1) - t(j);
    % % % % %             else
    % % % % %                 dt = t(j)-t(j-1);
    % % % % %             end
    % % %             ys(:,j) = expm(As(:,:,i)*dt)*ys(:,j-1) + As(:,:,i)^-1*(expm(As(:,:,i)*dt) - eye(n))*Bs(:,:,i)*u;
    % % % %     %         ys(1,j) = min(max(ys(1,j), -2), u(1)+2);
    % % %         end
        end
    end

    ys = [Xss(:,1) ys];
    
    if sum(isnan(ys))
        ME = MException('resultisNaN:noSuchVariable', ...
                           'Waveform reconstruction resulted in NaN');
        throw(ME);
    end
catch ERROR
    ERROR
end

end
