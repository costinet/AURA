function [JoutStart,JoutEnd] = discreteJacobianConstraint(obj)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    Cbnd = obj.topology.Cbnd(:,:,obj.converter.swind);

    [Jt, ~] = discreteJacobian(obj, 2);
        
    Jout = zeros(size(Cbnd,1),size(Jt,2),size(Jt,3));
    JoutStart = Jout;
    JoutEnd = Jout;
    %J(state, at time, time changed)
    for i = 1:size(Jt,3)
        for j = 1:size(Jt,2)
            Jout(:,j,i) = Cbnd(:,:,j)*Jt(:,j,i);

            % Cbnd is during the interval, but Jt is about the states at the
            % interface.
            JoutStart(:,j,i) = Cbnd(:,:,j)*Jt(:,j,i);
            endInt = circshift(1:size(Jt,2),-1);
            JoutEnd(:,j,i) = Cbnd(:,:,j)*Jt(:,endInt(j),i);
        end
    end
end

