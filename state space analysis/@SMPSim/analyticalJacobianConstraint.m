function [JoutStart,JoutEnd] = analyticalJacobianConstraint(obj)
%analyticalJacobianConstraint calculate discrete Jacobian with respect to
%constraints as output signals using the analytical jacobian
%   
%   See Also SMPSim.analyticalJacobian, SMPSim.discreteJacobianConstraint

    Cbnd = obj.topology.Cbnd(:,:,obj.converter.swind);

    [Jt, ~] = analyticalJacobian(obj);
    endInt = circshift(1:size(Jt,2),-1);
        
    JoutStart = zeros(size(Cbnd,1),size(Jt,2),size(Jt,3));
    JoutEnd = JoutStart;
    for i = 1:size(Jt,3)    
        %J(state, at time, time changed)
        %Cbnd(constraint, state, at time)
        JoutStart(:,:,i) = pagemtimes(Cbnd,permute(Jt(:,:,i),[1 3 2]));
        JoutEnd(:,:,i) = pagemtimes(Cbnd,permute(Jt(:,endInt,i),[1 3 2]));
    end
end
