function [J, J2, XssF, XssB, X0, dt] = discreteJacobian(obj, order)
%discreteJacobian calculate discrete Jacobian numerically
%   
%   [J] = discreteJacobian(obj)
%   for the SMPSim object obj, calculates discrete Jacobian J.  The
%   Jacobian is a 3-dimensional matrix with structure J(i, j, k): partial
%   derivative of X(i,j) due to perturbations in t(k). By default, the
%   function uses first-order forward difference approximation.  The time
%   perturnation dT is dictated by SMPSconverter.getDeltaT
%
%   [J, J2] = discreteJacobian(obj, order)
%   optional variable order can be either 1 or 2.  If set to two, a
%   second-order central differences approximation is used to generate J2.
%   J is unaffected.
%
%   [J, J2, XssF, XssB, X0, dt] = discreteJacobian(...)
%   optional outputs XssF, XssB, and X0 are the forward, backward and 
%   nominal state vectors at each time interval.  The calculated timestep
%   for the differences is dt.
%   
%
%   See Also SMPSconverter.getDeltaT, SMPSim.FindValidSteadyState

   
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
        dtApp = thats(i);

        tps = validateTimePerturbations(obj.converter, i, dtApp);
        XssF = obj.steadyState(tps);%perturbedSteadyState(obj, tps);
        dt1 = tps(i);
        dt = dt1;

        
        if order == 1
            J(:,:,i) = (XssF(:,1:end-1)-X0(:,1:end-1))/dt1;
            J2 = [];
            XssB = [];
        elseif order == 2

            tps = validateTimePerturbations(obj.converter, i, -dtApp);
            XssB = obj.steadyState(tps);%perturbedSteadyState(obj, tps);
            dt2 = tps(i);
            
            dt = [dt1, dt2];
            
            J(:,:,i) = (XssF(:,1:end-1)-XssB(:,1:end-1))/(dt1-dt2);
            J2(:,:,i) = (XssF(:,1:end-1) - 2*X0(:,1:end-1) + XssB(:,1:end-1))/((dt1-dt2)/2)^2;
        end
        
    end
    
    assert(all(abs(obj.ts - oldts) < obj.converter.timingThreshold), 'Unintended modification to timing during jacobian');

        
end

