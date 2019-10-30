function [fault,time_changed] = fast_dynamic_check(obj,i,j,k,Xs)
%fast_dynamic_check Takes the fastest natural frequency of the cirucit
%and check to see if there is a violation within the subinterval
%defined by i, j, k
%   Detailed explanation goes here

%     _   _   _  ____    _
%    / \ | | | |/ _  |  / \
%   / _ \| | | | (_| | / _ \
%  / ___ | |_| |> _  |/ ___ \
% /_/   \_\___//_/ |_/_/   \_\


fault = 0;
time_changed = [];
for n = 1:1:length(obj.eigA(:,j-1))
    fault = 0;
    if sign(imag(obj.eigA(n,j-1)))==1 && obj.ts(j-1)>2.*pi/imag(obj.eigA(n,j-1))
        
        % Value at 0
        X_quarters = Xs(:,j-1);
        xdot_quarters = obj.As(:,:,k)*X_quarters+obj.Bs(:,:,k)*obj.u;
        
        changed = [];
        
        
        for halfes = 1:1:4
            
            % Value of pi/2 pi 3pi/2 2pi
            M = expm(halfes*pi()/2/imag(obj.eigA(n,j-1))*[obj.As(:,:,k), eye(length(obj.As(:,:,k))); zeros(size(obj.As(:,:,k),1), 2*size(obj.As(:,:,k),2))]);
            intEAt = M(1:size(obj.As(:,:,k),1), size(obj.As(:,:,k),2)+1:end);
            fresp = intEAt*(obj.Bs(:,:,k)*obj.u);
            X_quarters(:,halfes+1)=expm(halfes*pi()/2/imag(obj.eigA(n,j-1))*obj.As(:,:,k))*Xs(:,j-1)+fresp;
            xdot_quarters(:,halfes+1) = obj.As(:,:,k)*X_quarters(:,halfes+1)+obj.Bs(:,:,k)*obj.u;
            
            if sign(xdot_quarters(i,halfes))==-1 && sign(xdot_quarters(i,halfes+1))==1
                changed = halfes;
            end
        end
        
        
        if isempty(changed)
            fault = 0;
            continue
        end
        
        left = (changed-1)*pi()/2;
        right = changed*pi()/2;
        
        
        while left<=right
            mid = (left+right)/2;
            
            M = expm(mid/imag(obj.eigA(n,j-1))*[obj.As(:,:,k), eye(length(obj.As(:,:,k))); zeros(size(obj.As(:,:,k),1), 2*size(obj.As(:,:,k),2))]);
            intEAt = M(1:size(obj.As(:,:,k),1), size(obj.As(:,:,k),2)+1:end);
            fresp = intEAt*(obj.Bs(:,:,k)*obj.u);
            X_mid=expm(mid/imag(obj.eigA(n,j-1))*obj.As(:,:,k))*Xs(:,j-1)+fresp;
            xdot_mid = obj.As(:,:,k)*X_mid+obj.Bs(:,:,k)*obj.u;
            
            if X_mid(i) <= -1
                fault  = 1;
                break
            end
            
            if sign(xdot_mid(i))==-1 && xdot_mid(i)*(2*pi()/imag(obj.eigA(n,j-1)) - mid/imag(obj.eigA(n,j-1)))+X_mid(i) > -1
                fault = 0;
                break
            end
            if sign(xdot_mid(i))==1
                right = mid;
            else
                left = mid;
                
            end
            
        end
    end
    if fault == 1
        
        Spacing = [0 pi()/2 pi() 3*pi()/2 2*pi()];
        Spacing = [Spacing(1:changed) mid Spacing(changed+1:end)];
        DT_points = [X_quarters(:,1:changed) X_mid X_quarters(:,changed+1:end)];
        index = find(DT_points(i,:) < -1, 1, 'first');
        
        delta_w = (-1-DT_points(i,index-1))/(DT_points(i,index)-DT_points(i,index-1))*(Spacing(i,index)-Spacing(i,index-1));
        
        step_needed = Spacing(i,index-1)+delta_w;
        
        time_changed = step_needed/2/pi()*2.*pi/imag(obj.eigA(n,j-1));
    end
end


end

