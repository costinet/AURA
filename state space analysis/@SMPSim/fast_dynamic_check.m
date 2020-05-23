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

try
    
    fault = 0;
    time_changed = [];
    
    Vf = 1;
    
    
    for n = 1:1:length(obj.eigA(:,j-1)) % Iterate through eigenvalues that are needed
        fault = 0;
        
        % This only advances if the eigen values are postive and the
        % dynamics of the time interval are greater than the one over
        % 2*frquency of the eigenvalues
        
        if sign(imag(obj.eigA(n,j-1)))==1 && obj.ts(j-1)>pi/2/imag(obj.eigA(n,j-1))
            
            % Value at 0
            X_quarters = Xs(:,j-1);
            xdot_quarters = obj.As(:,:,k)*X_quarters+obj.Bs(:,:,k)*obj.u;
            
            changed = [];
            
            if  obj.ts(j-1)>2*pi/imag(obj.eigA(n,j-1))
                iterations = 4;
            else
                iterations = floor(obj.ts(j-1)/(pi/2/imag(obj.eigA(n,j-1))));
            end
            
            for halfes = 1:1:iterations
                
                % Value and derivative of the time interval at pi/2 pi 3pi/2 2pi of the eigenvalue being observed
                M = expm(halfes*pi()/2/imag(obj.eigA(n,j-1))*[obj.As(:,:,k), eye(length(obj.As(:,:,k))); zeros(size(obj.As(:,:,k),1), 2*size(obj.As(:,:,k),2))]);
                intEAt = M(1:size(obj.As(:,:,k),1), size(obj.As(:,:,k),2)+1:end);
                fresp = intEAt*(obj.Bs(:,:,k)*obj.u);
                X_quarters(:,halfes+1)=expm(halfes*pi()/2/imag(obj.eigA(n,j-1))*obj.As(:,:,k))*Xs(:,j-1)+fresp;
                xdot_quarters(:,halfes+1) = obj.As(:,:,k)*X_quarters(:,halfes+1)+obj.Bs(:,:,k)*obj.u;
                
                % Checks to see if there was a sign change in derivative
                % of any of the state variables from negative to
                % positive to indicate that there could be a point
                % where the diode should have turned ON at some point
                % between the discrete time points
                if sign(xdot_quarters(i,halfes))==-1 && sign(xdot_quarters(i,halfes+1))==1
                    changed = halfes;
                end
            end
            
            % Builds the value and derivative at the fastest
            % eigenvalue frequency
            if iterations < 4
                X_quarters(:,halfes+2) = Xs(:,j);
                xdot_quarters(:,halfes+2) = obj.As(:,:,k)*X_quarters(:,halfes+2)+obj.Bs(:,:,k)*obj.u;
                
                if sign(xdot_quarters(i,halfes+1))==-1 && sign(xdot_quarters(i,halfes+2))==1
                    changed = halfes+1;
                end
            end
            % If there is no sign change then continue on with the
            % interations: they dynamics of this eigenvalue are slower
            % than the interval observed
            if isempty(changed)
                fault = 0;
                continue
            end
            
            % If there was an observed change in derivative of the state
            % variables then need to determine if there is a violation at
            % the local min or max in the time interval observed
            
            
            % Debug code to determine why there is a difference between
            % projection when their is a violation becoming positive
            
            
            %         debug =1;
            % if debug
            % T_fast = 2*pi/imag(obj.eigA(n,j-1));
            % SS = ss(obj.As(:,:,k), obj.Bs(:,:,k), eye(size(obj.As,1)), 0);
            %     tsim = linspace(0, T_fast, 1000);
            %
            %     if ts(Ti)
            %         [~, ~, x] = lsim(SS, obj.u*ones(size(tsim)), tsim, Xs(:,Xic-1));
            %         xdot = As(:,:,Xic-1)*x' + Bs(:,:,Xic-1)*u;
            %     else
            %         x = Xs(:,Ti)';
            %         xdot = As(:,:,Xic-1)*x' + Bs(:,:,Xic-1)*u;
            %     end
            %
            %
            % end
            %
            
            
            
            
            
            
            
            % If the last negative derivative point is positive
            if X_quarters(i,changed) > -Vf
                T_fast = 2*pi/imag(obj.eigA(n,j-1));
                % Set up binary search
                left = (changed-1)*pi()/2;
                right = changed*pi()/2;
                
                t_mid=(((-Vf)-X_quarters(i,changed))/xdot_quarters(i,changed))+T_fast*left/2/pi();
                if t_mid>T_fast
                    fault = 0;
                    continue
                end
                
                mid = t_mid*2*pi()/T_fast;
                
                while left<=right
                    
                    % This is not a true binary search where the middle value
                    % is always picked after every iteration, instead the left
                    % post value is projected to -1
                    
                    M = expm(mid/imag(obj.eigA(n,j-1))*[obj.As(:,:,k), eye(length(obj.As(:,:,k))); zeros(size(obj.As(:,:,k),1), 2*size(obj.As(:,:,k),2))]);
                    intEAt = M(1:size(obj.As(:,:,k),1), size(obj.As(:,:,k),2)+1:end);
                    fresp = intEAt*(obj.Bs(:,:,k)*obj.u);
                    X_mid=expm(mid/imag(obj.eigA(n,j-1))*obj.As(:,:,k))*Xs(:,j-1)+fresp;
                    xdot_mid = obj.As(:,:,k)*X_mid+obj.Bs(:,:,k)*obj.u;
                    
                    
                    % Change the 
                    if X_mid(i) <= -Vf
                        fault  = 1;
                        break
                    end
                    
                    if xdot_mid(i) >0
                        fault = 0;
                        break
                    end
                    
                    if sign(xdot_mid(i))==-1 && xdot_mid(i)*(2*pi()/imag(obj.eigA(n,j-1)) - mid/imag(obj.eigA(n,j-1)))+X_mid(i) >= -1
                        fault = 0;
                        break
                    end
                    if sign(xdot_mid(i))==1
                        right = mid;
                    else
                        left = mid;
                        
                    end
                    
                    if 0>(((-1)-X_mid(i))/xdot_mid(i))
                        J = 441564514;
                    end
                    
                    t_mid=(((-1)-X_mid(i))/xdot_mid(i))+T_fast*left/2/pi();
                    if t_mid>T_fast
                        fault = 0;
                        break
                    end
                    
                    if t_mid<0
                        J = 8578923045934;
                    end
                    

                    if (((-1)-X_mid(i))/xdot_mid(i)) < 1e-15 && sign(xdot_mid(i))==-1
                        fault = 1;
                        break
                    end
                    
                    
                    mid = t_mid*2*pi()/T_fast;
                end
                
                
                
                % If the last negative derivative point is negative (beyond the forward voltage)
            elseif X_quarters(i,changed) < -1
                changed  = changed - 1;
                if changed  ==0
                    mid = 0;
                    break
                end
                T_fast = 2*pi/imag(obj.eigA(n,j-1));
                % Set up binary search
                left = (changed-1)*pi()/2;
                right = changed*pi()/2;
                
                t_mid=(((-1)-X_quarters(i,changed))/xdot_quarters(i,changed))+T_fast*left/2/pi();
                if t_mid>T_fast
                    fault = 0;
                    continue
                end
                
                mid = t_mid*2*pi()/T_fast;
                
                while left<=right
                    
                    % This is not a true binary search where the middle value
                    % is always picked after every iteration, instead the left
                    % post value is projected to -1
                    
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
                    
                    t_mid=(((-1)-X_quarters(i,changed))/xdot_quarters(i,changed))+T_fast*left/2/pi();
                    if t_mid>T_fast
                        fault = 0;
                        continue
                    end
                    
                    mid = t_mid*2*pi()/T_fast;
                end
            end
            
            
        end
        
        if fault == 1
            
            % This is not needed any more since I am not doing a pure
            % binary search.
            
            %{
        % Places the radian value that is close to the transistion in
        % the Spacing array of radians from 0 to 2pi
        Spacing = [0 pi()/2 pi() 3*pi()/2 2*pi()];
        Spacing = [Spacing(1:changed) mid Spacing(changed+1:end)];
        
        DT_points = [X_quarters(:,1:changed) X_mid X_quarters(:,changed+1:end)];
        
        index = find(DT_points(i,2:end) < -1, 1, 'first');
        
        delta_w = (-1-DT_points(i,index-1))/(DT_points(i,index)-DT_points(i,index-1))*(Spacing(index)-Spacing(index-1));
        
        step_needed = Spacing(index-1)+delta_w;
        
        time_changed = step_needed/2/pi()*2.*pi/imag(obj.eigA(n,j-1));
            %}
            
            
            time_changed = mid/imag(obj.eigA(n,j-1));
            
        end
    end
    
catch ME
    disp(ME);
end
end

