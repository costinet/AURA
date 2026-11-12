function Err = compareJacobians(obj)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    obj SMPSim
end

    [Ja, dta] = analyticalJacobian(obj);
    % [J, J2, XssF, XssB, X0, dt] = discreteJacobian(obj, 1)

    Xs = obj.steadyState();
    Jd = zeros(size(Ja));
    for i = 1:length(obj.ts)
        dts = 0*obj.ts;
        dts(i) = min(obj.ts)/100000;
        dXs(:,:,i) = obj.steadyState(dts);

        for j = 1:size(dXs,2)
            dXs(:,j,i) = obj.Is(:,:,i)*dXs(:,j,i);
        end
        
        Jd(:,:,i) = (dXs(:,1:end-1,i)-Xs(:,1:end-1))./dts(i);
    end

    Err = norm(Jd-Ja,'fro')./Jd.^2;

end