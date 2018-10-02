function [fresp, intEAt] = forcedResponse(obj, A, expA, B, u, t, storeResult)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

    ns = size(B,1);
    
    if(nargin ==6)
        storeResult = 0;
    end
    
    % Check if this has already been calculated
    if(storeResult)
        index = find(sum(sum(bsxfun(@eq,obj.oldAs, A))) == numel(A),1,'first');
        if(index)
            index = index*(obj.oldts(index) == t);
        end
    else
        index = 0;
    end
    
    if index
        intEAt = obj.oldIntEAt(:,:,index);
    elseif(cond(full(A)) < obj.condThreshold) % not ill-conditioned, use matrix inverse
        intEAt = A\(expA-eye(ns));
    else % As(i) is ill-conditioned
       M = expm(t*[A, eye(length(A)); zeros(size(A,1), 2*size(A,2))]);
       intEAt = M(1:size(A,1), size(A,2)+1:end);
    end
    fresp = intEAt*(B*u);
end

