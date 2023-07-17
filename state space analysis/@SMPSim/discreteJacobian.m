function [J, J2, XssF, XssB, X0, dt] = discreteJacobian(obj, order)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%   J(x, y, z) = x-> which state, y-> at which time interval , z->which time interval perturbed
   
    assert(~isempty(obj.Xs), 'steady-state solution must be found before computing jacobian');
    
    nt = size(obj.Xs,2)-1;
    ns = size(obj.Xs,1);
    
    J = zeros(ns, nt, nt);
    
    
    if nargin == 1
        order = 1;
    else
        J2 = zeros(ns, nt, nt);
    end
    
    oldts = obj.ts; %store timing
    X0 = obj.Xs;    %store original steady-state
    
    [thats, ~] = getDeltaT(obj.converter);
    thats = min([thats; circshift(thats,-1)],[],1);
    
    for i = 1:nt % each switching interval
%         t0 = obj.ts(i);
%         if i<nt
%             t1 = obj.ts(i+1);
%         else
%             t1 = obj.ts(1);
%         end
%         
%         if i>1
%             t2 = obj.ts(i-1);
%         else
%             t2 = obj.ts(end);
%         end    
%         dt = min([sum(obj.ts)/1000, t0/2, t1/2, t2/2]);
        dtApp = thats(i);
        
        
%         newts = obj.ts;
%         newts(i) = newts(i) + dt;
%         if i<nt
%             newts(i+1) = newts(i+1) - dt;
%         else
%             newts(1) = newts(1) - dt;
%         end
%         obj.ts = newts;

% % %         dt = adjustUncontrolledTiming(obj.converter, i, dt);
% % %         XssF = obj.steadyState;
% % % %         adjustUncontrolledTiming(obj.converter, i, -dt);
% % %         obj.converter.undoLatestTimeChange

        tps = validateTimePerturbations(obj.converter, i, dtApp);
        XssF = obj.steadyState(tps);%perturbedSteadyState(obj, tps);
        dt1 = tps(i);
        dt = dt1;

        
        if order == 1
            J(:,:,i) = (XssF(:,1:end-1)-X0(:,1:end-1))/dt1;
            J2 = [];
            XssB = [];
        elseif order == 2
%             newts = oldts;
%             newts(i) = newts(i) - dt;
%             if i<nt
%                 newts(i+1) = newts(i+1) + dt;
%             else
%                 newts(1) = newts(1) + dt;
%             end
%             obj.ts = newts;
% % %             dt = adjustUncontrolledTiming(obj.converter, i, -dt);
% % %             XssB = obj.steadyState;
% % % %             adjustUncontrolledTiming(obj.converter, i, dt);
% % %             obj.converter.undoLatestTimeChange

            tps = validateTimePerturbations(obj.converter, i, -dtApp);
            XssB = obj.steadyState(tps);%perturbedSteadyState(obj, tps);
            dt2 = tps(i);
            
            dt = [dt1, dt2];
            
            J(:,:,i) = (XssF(:,1:end-1)-XssB(:,1:end-1))/(dt1-dt2);
            J2(:,:,i) = (XssF(:,1:end-1) - 2*X0(:,1:end-1) + XssB(:,1:end-1))/((dt1-dt2)/2)^2;
        end
        
%         obj.ts = oldts; %restore original timing
    end
    
    assert(all(abs(obj.ts - oldts) < obj.converter.timingThreshold), 'Unintended modification to timing during jacobian');

        
end

