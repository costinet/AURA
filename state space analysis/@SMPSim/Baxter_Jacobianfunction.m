function  [J, J2, XssF, XssB, X0, dt] = Baxter_Jacobianfunction(obj, order)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%   J(x, y, z) = x-> which state, y-> at which time interval , z->which time interval perturbed
   
    assert(~isempty(obj.Xs), 'steady-state solution must be found before computing jacobian');
    keep_SS = false;
    varToPerturb = 'ts';
    nt = size(obj.Xs,2)-1;
    ns = size(obj.Xs,1);
    
    J = zeros(ns, nt, nt);
    ts = obj.ts;
    
    if nargin == 1
        order = 1;
    else
        J2 = zeros(ns, nt, nt);
    end
    
    pI = 1:nt;
    cI = circshift(pI,1);
    X0 = obj.Xs;
    
    delta_DTs =10e-14; %max(min(ts)/1000, sum(ts)/100000); 
            if  max(eps(ts))*1000>delta_DTs ||  min(pi./1000./abs(imag(obj.eigA(:,:))),[],'all') < delta_DTs
                fprintf('not good')
                
            end
    
    for i = 1:nt % each switching interval
       [XssF,dt1] = Baxter_StateSensitivity2(obj,keep_SS, varToPerturb, pI(i), delta_DTs, cI(i));
       
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
            XssB = perturbedSteadyState(obj, tps);
            dt2 = tps(i);
            
            dt = [dt1, dt2];
            
            J(:,:,i) = (XssF(:,1:end-1)-XssB(:,1:end-1))/(dt1-dt2);
            J2(:,:,i) = (XssF(:,1:end-1) - 2*X0(:,1:end-1) + XssB(:,1:end-1))/((dt1-dt2)/2)^2;
        end
        
%         obj.ts = oldts; %restore original timing
    end
    
    %assert(all(abs(obj.ts - oldts) < obj.converter.timingThreshold), 'Unintended modification to timing during jacobian');

        
end

