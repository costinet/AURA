function [J, J2, XssF, XssB, X0, dt] = discreteJacobian(obj, order)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
   
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
    
    for i = 1:nt % each switching interval
        t0 = obj.ts(i);
        if i<nt
            t1 = obj.ts(i+1);
        else
            t1 = obj.ts(1);
        end
        dt = min([sum(obj.ts)/1000, t0/2, t1/2]);
        newts = obj.ts;
        newts(i) = newts(i) + dt;
        if i<nt
            newts(i+1) = newts(i+1) - dt;
        else
            newts(1) = newts(1) - dt;
        end
        
        obj.ts = newts;
        XssF = obj.SS_Soln();
        
        if order == 1
            J(:,:,i) = (XssF(:,1:end-1)-X0(:,1:end-1))/dt;
        elseif order == 2
            newts = oldts;
            newts(i) = newts(i) - dt;
            if i<nt
                newts(i+1) = newts(i+1) + dt;
            else
                newts(1) = newts(1) + dt;
            end
            
            obj.ts = newts;
            XssB = obj.SS_Soln();
            
            J(:,:,i) = (XssF(:,1:end-1)-XssB(:,1:end-1))/2/dt;
            J2(:,:,i) = (XssF(:,1:end-1) - 2*X0(:,1:end-1) + XssB(:,1:end-1))/dt^2;
        end
        
        obj.ts = oldts; %restore original timing
    end
    
    
        
        
end

