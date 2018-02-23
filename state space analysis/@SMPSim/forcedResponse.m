function [fresp, intEAt] = forcedResponse(obj, A, expA, B, u, t, storeResult)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

    ns = size(B,1);
    
    if(nargin ==6)
        storeResult = 0;
    end
    
    % Check if this has already been calculated
    if(storeResult)
        index = find(sum(sum(obj.oldAs == A)) == numel(A),1,'first');
    else
        index = 0;
    end
    
    if index
        intEAt = obj.oldIntEAt(:,:,index);
    elseif(cond(A) < obj.condThreshold) % not ill-conditioned, use matrix inverse
        intEAt = A^-1*(expA-eye(ns));
    else % As(i) is ill-conditioned
       intEAt = zeros(ns,ns);

       % Try iterative summation
       % NOTE: may be useful to store Ai^k for 1<k<100 in class-based
       % implementation
       k = 0;
       while sum(sum(~isfinite(A^k))) == 0 && k < 100
           convError = 1/prod(1:k+1)*A^k*t^(k+1);
           intEAt = intEAt + convError;           
           k = k+1;
       end

        % Didn't work, try lsim approach (numerical integration)
        if sum(sum(isnan(intEAt))) || sum(sum(isinf(intEAt))) || max(max(convError))>1e-9
           x0 = zeros(1,ns)';
           for i=1:ns
               Btest = zeros(ns,1);
               Btest(i) = 1;
               [~,~,test3] = lsim(ss(A,Btest,zeros(1,ns),0),ones(1,100), linspace(0,t,100),x0);
               intEAt(:,i) = test3(end,:)';
           end

           % If integration had problems, try symbolic Laplace (very slow; last resort)
           if sum(sum(isnan(intEAt))) || sum(sum(isinf(intEAt)))
               syms s t
               intEAt = eval(subs(vpa(ilaplace(((s*eye(ns)-A)^-1)/s)), t, t));
           end
        end  
    end
    fresp = intEAt*(B*u);
end

