function [ avgXs, avgYs ] = ssAvgs(obj, Xss)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here

    As = obj.As;
    Bs = obj.Bs;
    Cs = obj.Cs;
    Ds = obj.Ds;
    ts = obj.ts;
    u = obj.u;

    invertible = 1;

    n = size(As,3);
    ns = size(Bs,1);
    nc = size(Cs,1);
    
    avgxs = zeros(ns,n);
    avgys = zeros(nc,n);

    for i=1:n
        if cond(As(:,:,i)) > obj.condThreshold
            invertible = 0;
            break
        end
    end

    if invertible
        for i=1:n
            avgxs(:,i) = inv(As(:,:,i))/ts(i) * ( (expm(As(:,:,i)*ts(i))-eye(ns))*Xss(:,i) + ...
                (As(:,:,i)\(expm(As(:,:,i)*ts(i))-eye(ns))-ts(i)*eye(ns))*Bs(:,:,i)*u);
            avgys(:,i) = Cs(:,:,i)*avgxs(:,i) + Ds(:,:,i)*u;
        end
        avgXs = sum(avgxs.*repmat(ts,ns,1),2)/sum(ts);
        avgYs = sum(avgys.*repmat(ts,nc,1),2)/sum(ts);  
    else 
        % double number of states
        As(2*ns, 2*ns, 1) = 0; 
        Bs(2*ns,1,1) = 0;
        Xss(2*ns,1,1) = 0;
        
        % Use new states to integrate exisitng states
        for i = 1:ns
            As(ns+i, i, :) = ones(1,n);
        end

        for i = 1:n
            Asp = sparse(As(:,:,i));
            [~, intEAt] = obj.forcedResponse(Asp, expm(Asp*ts(i)), Bs(:,:,i), u, ts(i));
            intStates = expm(Asp*ts(i))*Xss(:,i) + intEAt*Bs(:,:,i)*u;
            avgxs(:,i) = intStates(ns+1:end)/ts(i);
            avgys(:,i) = Cs(:,:,i)*avgxs(:,i) + Ds(:,:,i)*u;
        end
        avgXs = sum(avgxs.*repmat(ts,ns,1),2)/sum(ts);
        avgYs = sum(avgys.*repmat(ts,nc,1),2)/sum(ts);  
    end


end

