function [outputArg1,outputArg2] = FindZeroDeriv(obj,i,j,k,debug)
%FINDZERODERIV Searches to find there is a zero derivative of the passed
%measurement in the give time interval
%   Detailed explanation goes here

%     %%%%%%   %      %  %%%%%%%    %%%%%%
%    %      %  %      %  %      %  %      %
%    %      %  %      %  %      %  %      %
%    %%%%%%%%  %      %  %%%%%%%   %%%%%%%%
%    %      %  %      %  %%        %      %
%    %      %  %      %  % %       %      %
%    %      %  %      %  %  %      %      %
%    %      %  %      %  %   %     %      %
%    %      %   %    %   %    %    %      %
%    %      %    %%%%    %     %   %      %


% Set variables from class
Xs = obj.Xs;
ts = obj.ts;

% Set varaibles for fmincon
I_solemnly_swear_that_I_am_up_to_no_good = [i,j,k];
A = [];
B = [];
Aeq = [];
Beq = [];
LB = 0;
UB = ts(j);
fun = @(t) t;

% Fix options for fmincon
options = optimoptions('fmincon','StepTolerance',1e-12,'ConstraintTolerance',1e-9);
if debug
    options = optimoptions('fmincon','Display','off','StepTolerance',1e-65,'ConstraintTolerance',1e-3,'OptimalityTolerance',1e-5);
end

% Run 1st simulation with x0 close to beginning of time interval
x0 = 0.1*ts(j);
[X1,FVAL1,EXITFLAG1,OUTPUT1,LAMBDA1,GRAD1,HESSIAN1] = fmincon(fun,x0,A,B,Aeq,Beq,LB,UB,@obj.zero,options);


% Run 2nd simulation with x0 close to end of time interval
x0 = 0.9*ts(j);
[X2,FVAL2,EXITFLAG2,OUTPUT2,LAMBDA2,GRAD2,HESSIAN2] = fmincon(fun,x0,A,B,Aeq,Beq,LB,UB,@obj.zero,options);

if (EXITFLAG1==1) && (EXITFLAG2==1) % If a solution is found for both fmincons

    if abs(X2-X1)>(ts(j)*0.5) % If there are multiple zero crossing (determined by numerical approximation that the same zero crossing will not numerically approximate to greate than 5% of time interval)
        
    
    
if debug
    J = zeros(size(Xs,1),size(Xs,1));
    J(i,i) = 1;
    syms t
    stack = [];
    t1 = linspace(0,5e-6,1000);
    for counts = 1:1:length(t1)
        eqn2 = J*(As(:,:,k)*expm(As(:,:,k).*t1(counts))*Xs(:,j)+expm(As(:,:,k).*t1(counts))*Bs(:,:,k)*u);
        stack(end+1) = eqn2(1,:);
    end
    plot(t1,stack)
end

end % That's all Folks

