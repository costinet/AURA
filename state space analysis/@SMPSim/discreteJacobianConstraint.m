function [JoutStart,JoutEnd] = discreteJacobianConstraint(obj)
%discreteJacobianConstraint calculate discrete Jacobian with respect to
%constraints as output signals
%   
%   See Also SMPSim.discreteJacobian

    Cbnd = obj.topology.Cbnd(:,:,obj.converter.swind);

    [Jt, ~] = discreteJacobian(obj, 2);
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

%% Implementation without pagemtimes
%     Cbnd = obj.topology.Cbnd(:,:,obj.converter.swind);
% 
%     [Jt, ~] = discreteJacobian(obj, 2);
%         
%     Jout = zeros(size(Cbnd,1),size(Jt,2),size(Jt,3));
%     JoutStart = Jout;
%     JoutEnd = Jout;
%     %J(state, at time, time changed)
%     for i = 1:size(Jt,3)
%         for j = 1:size(Jt,2)
%             Jout(:,j,i) = Cbnd(:,:,j)*Jt(:,j,i);
% 
%             % Cbnd is during the interval, but Jt is about the states at the
%             % interface.
%             JoutStart(:,j,i) = Cbnd(:,:,j)*Jt(:,j,i);
%             endInt = circshift(1:size(Jt,2),-1);
%             JoutEnd(:,j,i) = Cbnd(:,:,j)*Jt(:,endInt(j),i);
%         end
%     end
% end

