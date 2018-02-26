function [ Xs] = SS_Soln( As, Bs, ts, u, Xi, Bi)
%HCM is the half-cycle matrix which negates any states which are odd with
%respect to the switching period.  For full cycle, just use HCM = eye(n)

if(nargin == 0)
    Dickson_solution_old
end
n = size(As,3);
ns = size(Bs,1);

condAs= zeros(1,n);

expAs = zeros(ns,ns,n);%[];
for i=1:n
    expAs(:,:,i) = expm(As(:,:,i)*ts(i));
    condAs(i) = cond(As(:,:,i));
end

cumProdExp = zeros(ns,ns,n);%[];
for i=1:n
    cumProdExp(:,:,i) = eye(ns);
    for k = 1:i
        cumProdExp(:,:,i) = expAs(:,:,k)*cumProdExp(:,:,i);
    end
end

cumProdExpRev = zeros(ns,ns,n);%[];
for i=1:n
    cumProdExpRev(:,:,i) = eye(ns);
    if(i<n)
        for k = i+1:n
            cumProdExpRev(:,:,i) = expAs(:,:,k)*cumProdExpRev(:,:,i);
        end
    end
end

RHSsum = zeros(ns,1);
for i=1:n
   if(condAs(i) < 1e9)
        RHSsum = RHSsum +  cumProdExpRev(:,:,i)*As(:,:,i)^-1*(expAs(:,:,i)-eye(ns))*(Bs(:,:,i)*u);
   else
       trickyTerm = zeros(ns,ns);
       k = 0;
       while sum(sum(~isfinite(As(:,:,i)^k))) == 0 && k < 100
           trickyTerm = trickyTerm + 1/factorial(k+1)*As(:,:,i)^k*ts(i)^(k+1);
           k = k+1;
       end
%        syms s t
%        test = eval(subs(ilaplace(((s*eye(ns)-As(:,:,i))^-1)/s), t, ts(i)));
%        test2 = expm(As(:,:,i)*ts(i))*eval(subs(ilaplace(((s*eye(ns)-As(:,:,i))^-1)/s), t, ts(i)))
%        trickyTerm - test
%        trickyTerm = test;
       RHSsum = RHSsum +  cumProdExpRev(:,:,i)*trickyTerm*(Bs(:,:,i)*u); 
   end
end


Xss = (eye(ns) - cumProdExp(:,:,n))^-1*RHSsum;

% % %% Hacky solution when result is off because (eye(ns) - cumProdExp(:,:,n)) is non-invertable
% % if(abs(Xss(1) - 12) > 12)
% %     warning(['Unable to find SS solution directly for Xi = ' num2str(Xi)]);
% %     Xss = Xi;%(eye(ns) - cumProdExp(:,:,n-1))^-1*RHSsum
% %     options = optimoptions('lsqlin','algorithm','trust-region-reflective','Display','none');
% %     A = (eye(ns) - cumProdExp(:,:,n));
% %     b = RHSsum;
% % %     w = diag(sqrt([1/12 1/24 1/36 1/5 1/9]));
% % %     A = w*A;
% % %     b = w*b;
% %     Xss = lsqlin(A, b, [],[],[],[], Xi.*Bi, Xi.*(2-Bi), Xss, options);
% % 
% % %     for i=1:50
% % %         err = (eye(ns) - cumProdExp(:,:,n))*Xss - RHSsum
% % % %         Xss = lsqlin((eye(ns) - cumProdExp(:,:,n)), RHSsum, [],[],[],[], Xss*.8, Xss*1.2)
% % %         Xss = Xss - [1 1 1 1 1]'.*err./Xss;
% % %     end
% % end

Xs(:,1) = Xss;
for(i=1:n)
   if(condAs(i) < 1e9)
        Xs(:,i+1) = expm(As(:,:,i)*ts(i))*Xs(:,i) + As(:,:,i)^-1*(expAs(:,:,i)-eye(ns))*(Bs(:,:,i)*u);
   else
       trickyTerm = zeros(ns,ns);
       k = 0;
       while sum(sum(~isfinite(As(:,:,i)^k))) == 0 && k < 100
           trickyTerm = trickyTerm + 1/factorial(k+1)*As(:,:,i)^k*ts(i)^(k+1);
           k = k+1;
       end
       Xs(:,i+1) = expm(As(:,:,i)*ts(i))*Xs(:,i) + trickyTerm*(Bs(:,:,i)*u); 
   end
end

if sum(isnan(Xss))
    ME = MException('resultisNaN:noSuchVariable', ...
                       'Waveform reconstruction resulted in NaN');
%     throw(ME);
end





