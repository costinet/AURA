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
       intEAt = zeros(ns,ns);
     
       %So far, this has never worked.
       % try diagonal modification of inverse
       Ap = A + eye(ns)*obj.gmin;
       k = 0;
       while(cond(Ap) > obj.condThreshold && k<10)
           Ap = A + eye(ns)*obj.gmin*2^k;
           k = k+1;
       end
       
       if(k<10)
           intEAt = Ap\(expA-eye(ns));
       else

           % Try iterative summation
           % NOTE: may be useful to store Ai^k for 1<k<100 
           % NOTE: significant convergence issues may need to be addressed.
           %  This seems to work if t^k/k! decreases faster than A^k blows
           %  up for k smaller than that where A^k overflows floating point
           %  values
           k = 0;
           while sum(sum(~isfinite(A^k))) == 0 && k < 100
               convError = 1/prod(1:k+1)*A^k*t^(k+1);
               intEAt = intEAt + convError;           
               k = k+1;
           end 
%            max(max(A))^(k-1)*t^k/factorial(k) <1  % estimate if it converged

            % Didn't work, try lsim approach (numerical integration)
            if sum(sum(isnan(intEAt))) || sum(sum(isinf(intEAt))) || max(max(convError))>1e-9
               fastDynamics = 1/(max(abs(eigs(A)))/2/pi);
               nsteps = min(t/fastDynamics*10,1e2);
               x0 = zeros(1,ns)';
               for i=1:ns
                   Btest = zeros(ns,1);
                   Btest(i) = 1;
                   [~,~,test3] = lsim(ss(full(A),Btest,zeros(1,ns),0),ones(1,nsteps), linspace(0,t,nsteps),x0);
                   intEAt(:,i) = test3(end,:)';
               end

               % If integration had problems, try symbolic Laplace (very slow; last resort)
               if sum(sum(isnan(intEAt))) || sum(sum(isinf(intEAt)))
                   ti = t;
                   syms s t
                   intEAt = eval(subs(vpa(ilaplace(((s*eye(ns)-A)^-1)/s)), t, ti));
               end
            end  
       end
    end
    fresp = intEAt*(B*u);
end

