function [ Xs] = SS_Soln(obj, Xi, Bi)
% Steady-state solution of switched system using state-space matrices
%
% [ Xs] = SS_Soln( As, Bs, ts, u, Xi, Bi) finds the state values Xs in
% steady-state for the switched system described by As, Bs, and ts
% according to
%
% dx/dt = Ai*x(t) + Bi*u ,
%
% for the ith interval.
%
% As is a 3-dimensional matrix of values for Ai, where As(:,:,i) is the 2D
% square matrix Ai during the ith interval.  
% Bs is a 3-dimensional matrix of values for Bi, where Bs(:,:,i) is the 2D
% matrix/vector Bi during the ith interval.
% ts is a vector of the time durations of each inverval
% u is the (assumed constant) independent input vector
% Xi is a vector "guess" of the results, which is used to determine
% solution validity in the case of singular matrices
% Bi is a vector of boundary coefficiencts used when iteration is required
% due to non-convergence of a solution.  A valid solution is required to be
% bounded within Bi.*Xi <= Xs <= (2-Bi).*Xi
%
% The result, Xs is a 2D matrix with ns rows and n+1 columns, where ns is the
% number of states in the system and n is the number of switching
% intervals.  The first and last column of Xs should be identical,
% corresponding to a valid steady-state solution.


As = obj.As;
Bs = obj.Bs;
ts = obj.ts;
u = obj.u;

if(nargin == 1)
    tryOpt = 0;
else
    tryOpt = obj.tryOpt;
end

n = size(As,3);
ns = size(Bs,1);

condAs= zeros(1,n);

expAs = zeros(ns,ns,n);
for i=1:n
    condAs(i) = cond(As(:,:,i));
    expAs(:,:,i) = expm(As(:,:,i)*ts(i));
end

cumProdExp = zeros(ns,ns,n);
for i=1:n
    cumProdExp(:,:,i) = eye(ns);
    for k = 1:i
        cumProdExp(:,:,i) = expAs(:,:,k)*cumProdExp(:,:,i);
    end
end

cumProdExpRev = zeros(ns,ns,n);
for i=1:n
    cumProdExpRev(:,:,i) = eye(ns);
    if(i<n)
        for k = i+1:n
            cumProdExpRev(:,:,i) = expAs(:,:,k)*cumProdExpRev(:,:,i);
        end
    end
end

RHSsum = zeros(ns,1);
trickyTerm = zeros(ns,1,n); %% A^-1(expm(A*T)-eye)*B*u
for i=1:n
   if(condAs(i) < 1e9) % not ill-conditioned, use matrix inverse
       trickyTerm(:,:,i) = As(:,:,i)^-1*(expAs(:,:,i)-eye(ns))*(Bs(:,:,i)*u);
%        RHSsum = RHSsum +  cumProdExpRev(:,:,i)*As(:,:,i)^-1*(expAs(:,:,i)-eye(ns))*(Bs(:,:,i)*u);
   else % As(i) is ill-conditioned
       iter = zeros(ns,ns);
       k = 0;
       
       % Try iterative summation
       % NOTE: may be useful to store Ai^k for 1<k<100 in class-based
       % implementation
       while sum(sum(~isfinite(As(:,:,i)^k))) == 0 && k < 100
%          convError = 1/factorial(k+1)*As(:,:,i)^k*ts(i)^(k+1); 
%          factorial appears to be slower than prod
           convError = 1/prod(1:k+1)*As(:,:,i)^k*ts(i)^(k+1);
           iter = iter + convError;           
           k = k+1;
       end
       trickyTerm(:,:,i) = iter*(Bs(:,:,i)*u);
       
       % Didn't work, try lsim approach (numerical integration)

       if sum(sum(isnan(trickyTerm(:,:,i)))) || sum(sum(isinf(trickyTerm(:,:,i)))) || max(max(convError))>1e-9
           x0 = ones(1,ns)';
           [~,~,test3] = lsim(ss(As(:,:,i),Bs(:,:,i),zeros(1,ns),0),repmat(u',100,1), linspace(0,ts(i),100),x0);
           trickyTerm(:,:,i) = test3(end,:)' - expAs(:,:,i)*x0;
           
           % If iterative summation didn't converge, try symbolic Laplace (very slow; last resort)
           if sum(sum(isnan(trickyTerm(:,:,i)))) || sum(sum(isinf(trickyTerm(:,:,i))))
               syms s t
               trickyTerm(:,:,i) = eval(subs(vpa(ilaplace(((s*eye(ns)-As(:,:,i))^-1)/s)), t, ts(i)))*(Bs(:,:,i)*u);
           end
       end  
   end
   RHSsum = RHSsum +  cumProdExpRev(:,:,i)*trickyTerm(:,:,i);
end


Xss = (eye(ns) - cumProdExp(:,:,n))^-1*RHSsum;

%% Hacky solution when result is off because (eye(ns) - cumProdExp(:,:,n)) is non-invertable
% Use optimization to solve without inversion through error minimization
% --> Still needs tweaking
if(tryOpt)
    if(abs(Xss(1) - Xi(1)) > 2)
        warning(['Unable to find SS solution directly for Xi = ' num2str(Xi)]);
        Xss = Xi;
        options = optimoptions('lsqlin','algorithm','trust-region-reflective','Display','none');
        A = (eye(ns) - cumProdExp(:,:,n));
        b = RHSsum;
        Xss = lsqlin(A, b, [],[],[],[], Xi.*Bi, Xi.*(2-Bi), Xss, options);
    end
end

%% from steady-state solution, go through and find states at each subinterval
Xs(:,1) = Xss;
for i=1:n
    Xs(:,i+1) = expAs(:,:,i)*Xs(:,i) + trickyTerm(:,:,i);
end

if sum(isnan(Xss))
    ME = MException('resultisNaN:noSuchVariable', ...
                       'Waveform steady-state solution resulted in NaN');
    throw(ME);
end

obj.Xs = Xs;


end

