function [ valid, newt, dist ] = validByInterval(obj, si, Xs)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    cnstr = obj.converter.topology.constraints;
    
    debug = 0;


    ns = length(obj.As(:,:,si));
    
    x0 = Xs(:,si);
    ti = obj.ts(si);
    
%     funt = @(t) -t;
    
%     fmincon(funt, ti, 

    
    
    
    %     [ x ] = stateValue_at_t( x0, t, si)
%             maxcheck(:,i) = cnstr.Cmax*Xs(:,i) + cnstr.Dmax*u - cnstr.Maxlim;
%         mincheck(:,i) = cnstr.Cmin*Xs(:,i) + cnstr.Dmin*u + cnstr.Minlim;
%         regcheck(:,i) = cnstr.Creg*Xs(:,i) + cnstr.Dreg*u - cnstr.regtarget;

    % find if state has gone above max limit
    funmax = @(t) [1 0]*((cnstr.Cmax*obj.stateValue_at_t( x0, t, si) + cnstr.Dmax*obj.u - cnstr.Maxlim).^2);
    
    % find if state has gone below min limit
    funmin = @(t) [1 0]*((cnstr.Cmin*obj.stateValue_at_t( x0, t, si) + cnstr.Dmin*obj.u - cnstr.Minlim).^2);
    
    % find zero-derivatives
    funderiv = @(t) ([1 0 0 0]*(obj.converter.topology.As(:,:,si)*obj.stateValue_at_t( x0, t, si) + ...
        obj.Bs(:,:,si)*obj.u)).^2;
    
    if(debug)
        figure(7); hold on;
        subplot(2,1,2); hold on;
        t = linspace(0, ti, 250);
        plot(t, funmax(t));
        plot(t, funmin(t));
        plot(t, funderiv(t)/10^16);
        
        subplot(2,1,1); hold on;
        plot(t, obj.stateValue_at_t( x0, t, si));
    end
    
    
    options = optimset('Display','off', 'Tolx', ti/1000); 
%         Tr = zeros(1, ns);
%     Tr(si) = 1;
%     A = obj.As(Ai);
%     funmax = @(x) (Tr*expm(A*x)*x0 - ub)^2;
%     funmin = @(x) (Tr*expm(A*x)*x0 - lb)^2;
    
%     [fresp, ~] = forcedResponse(obj, A, expA, B, u, t, storeResult)
    
    [X,F0,EXITFLAG,out]  = fminbnd(funderiv,0, ti,options);
    
    if(EXITFLAG == 1) % Found a zero-derivative
        if X > ti/100 && X<ti*99/100
            valid = 0;
            newt = X;
            dist = F0;
            return
        end
    else
        valid = 1;
        newt = -1;
    end

    [X,F1,EXITFLAG,out]  = fminbnd(funmax,0, ti,options);
    
    dist(1) = F1;
    
    %% If search converged to a point that isn't the end of the interval, or if the final value is over the max
    if(EXITFLAG == 1) &&  ((X<ti*99/100) || sum((cnstr.Cmax*Xs(:,si) + cnstr.Dmax*obj.u - cnstr.Maxlim)>0))
        valid = 0;
        newt = X;
    else
        valid = 1;
        newt = -1;
    end
    
%     [X,F2,EXITFLAG,~]  = fminbnd(funmin,0, ti,options);
%     dist(2) = F2;
%     
%      %% If search converged to a point that isn't the end of the interval, or if the final value is under the min
%     if(EXITFLAG == 1) &&  ((X<ti*99/100) || sum(cnstr.Cmin*Xs(:,si) + cnstr.Dmin*obj.u + cnstr.Minlim))
% %     if(EXITFLAG == 1)
%         valid = 0;
%         if(newt > 0)
%             newt = min(X,newt);
%         else
%             newt = X;
%         end
%     end
    
end

