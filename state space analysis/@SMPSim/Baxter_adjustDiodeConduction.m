function [ ts, dxsdt, hardSwNecessary, multcross, overresonant,keep_SS] = Baxter_adjustDiodeConduction(obj, Xs, Xic, Sir, Vmax, Vmin, progBar,Vio,keep_SS)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Xic = index to states (Xs) for discrete time interval which (may) contain
% errors due to diode conduction (column index of Xs)
% Sir = index to states (Xs) for the state in quesiton (row index of Xs)
% Vmax is the voltage maximum for the state
% Vmin is the voltage minimum for the state
% This changes depending on if it is a diode or a FET!!!!!!!!!!!!

As = obj.As;
Bs = obj.Bs;
ts = obj.ts;
u = obj.u;

%keep_SS = false;
change = [];
tsmax = [.7 .3 .3 .7]; % JUST FOR TESTING

try
    debug = 1;
    if(debug == 1)
        origts = ts;
    end
    
    if nargin < 7
        progBar = 0.01;
    end
    
    imagine_eigA = imag(obj.eigA); % get the imaginary part of the eigenvalue
    imagine_eigA(imagine_eigA(:)<=0) = 0; % Only take the postive part of the eigenvalue
    imagin_eigA_reduced = imagine_eigA(:,Xic-1:Xic);
    
    %maxStep =  min(min(pi()/5./imagine_eigA));
    
    maxStep =  min(min(pi()/5./imagin_eigA_reduced));
    % sum(ts)/100/(.5+1.5*progBar)
    
    if ~(Xic>1) || ~(Xic<=size(Xs,2)) % Check to see if Xic is valid
        ME = MException('resultisNaN:noSuchVariable', ...
            'Xi must be between 2 and length(Xs)');
        throw(ME);
    end
    
    
    Ti = Xic-1; % Current index of time interval
    Tc = Ti+1; % Next time interval
    
    % If the next time interval reaches the end of the defined period
    % then restart at beginning (1)
    if(Tc == length(ts) + 1)
        Tc = 1;
    end
    
    dxsdt = zeros(size(Xs));
    hardSwNecessary = 0;
    multcross = 0;
    overresonant = 0;
    %{
    % Run lsim for given time period
    C = zeros(1,size(As,1), 1);
    SS = ss(As(:,:,Xic-1), Bs(:,:,Xic-1), C, 0);
    tsim = linspace(0, ts(Ti), 100);
    
    if ts(Ti)
        [~, ~, x] = lsim(SS, u*ones(size(tsim)), tsim, Xs(:,Xic-1));
        xdot = As(:,:,Xic-1)*x' + Bs(:,:,Xic-1)*u;
    else
        x = Xs(:,Ti)';
        xdot = As(:,:,Xic-1)*x' + Bs(:,:,Xic-1)*u;
    end
    
    % Determine if there is a slope change in the state variable
    % and function being examined
    slopechange = find(diff(sign(xdot(Sir,:))));
    if ~isempty(slopechange)
        % This is a crude way to account for there being a little
        % bit of transient at the beginnining of a wavefrom that
        % is otherwise linear
        if abs(sum(sign(xdot(Sir,:))))>length(xdot(Sir,:))*0.95
            slopechange  = [];
        end
    end
    
    hardSwNecessary = 0;
    
    % adjust limits to account for possible I*Ron outside of bounds.
    %Vmax = max(Vmax, Xs(Sir, Xic-1));
    %Vmin = min(Vmin, Xs(Sir, Xic-1));
    
    % Find if the value of voltage exceeds the
    Vmaxcross = find(diff(x(:,Sir)>Vmax));
    Vmincross = find(diff(x(:,Sir)<Vmin));
    
    multcross = length(Vmaxcross) + length(Vmincross) > 1;
    overresonant = ~isempty(slopechange);
    massiveOvershoot = Xs(Sir,Xic) > Vmax + max(abs([Vmax Vmin 1])) || Xs(Sir,Xic) < Vmin - max(abs([Vmax Vmin 1]));
    
%     %if multcross || overresonant || massiveOvershoot
%         %% case where deadtime ringing causes multiple diode conduction intervals
%         % dead time is way too long, so set it to first crossing (unintelligently)
%         % and try again on next iteration
%         tdelta = tsim(max(min([Vmaxcross' Vmincross' slopechange]),2)) - ts(Xic-1);
%
%         if isempty(tdelta)
%             return
%         end
%
%         %guaranteed to be > 0, so no need to check
%         ts(Ti) = ts(Ti) + tdelta;
%         ts(Tc) = ts(Tc) - tdelta;
%         if(debug), disp('multcross || overresonant || massiveOvershoot found'); end
%         return;
%    % else
    %}
    %     if isempty(Vmaxcross) && isempty(Vmincross) && isempty(slopechange)
    %         % no diode conduction problem
    
    %% check for ZVS where possible
    %{
        There are some possible issues with the sensitivity beign too
        large when it is close to the zero crossing of the inductor.
        The sensitivity is too big, currently look to set sensitivity
        to be very small and see if resutls are okay-ish
        
        Test out using eps() and see if I run into any floating point
        issues
        
    %}
    
    delta_DTs =10e-14; %max(min(ts)/1000, sum(ts)/100000); % This is the area where there is an issue
    if  max(eps(ts))*1000>delta_DTs ||  min(pi./1000./abs(imag(obj.eigA(:,Xic)))) < delta_DTs
        fprintf('not good')
        
    end
    
    
    
    
    [dXs,delta_DTs] = obj.Baxter_StateSensitivity2(keep_SS, 'ts', Ti, delta_DTs, Tc);
    dxsdt = (dXs-Xs)/delta_DTs;
    [dXs_back,delta_DTs_back] = obj.Baxter_StateSensitivity2(keep_SS, 'ts', Tc, delta_DTs, Ti);
    dxsdt_back = (Xs - dXs_back)/delta_DTs_back;
    
    if delta_DTs_back == delta_DTs
        Central_Difference = (0.5.*dXs-0.5.*dXs_back)/(delta_DTs);
        if sum(sum((sign(dxsdt) == sign(dxsdt_back))-1)) ~= 0 || sum(sum((sign(dxsdt_back) == sign(Central_Difference))-1))~=0
            J = 1234123;
        end
    end
    
    if sum(sum((sign(dxsdt) == sign(dxsdt_back))-1)) ~= 0
        J = 54654853;
    end
    
    if   sign(dxsdt(Sir,Xic))~= sign(dxsdt_back(Sir,Xic))
        J = 456456;
    end
    
    Second_Derivative = (dXs_back-2*obj.Xs+dXs)/(delta_DTs^2);
    
    
    if(Vio==1)
        tdelta = (Vmax-Xs(Sir,Xic))/(dxsdt(Sir,Xic));
        %tdelta = min(max(tdelta, -ts(Ti) + delta_DTs), tsmax(Ti) - ts(Ti));
        tdelta = sign(tdelta)*min(abs(tdelta), maxStep);
       
        
        %  if(Xs(Sir,Xic)-Vmax < 0 && tdelta < 0)
        % Slope is positive, and non-ZVS, but partial tells me to reduce dead time
        % effect will be small (reduced current due to reduced t1)
        % so tdelta will always saturate and oscillate. Reduce step size
        % significantly to prevent oscillating
        %  tdelta = sign(tdelta)*min(abs(tdelta), maxStep/10);
        % end
        ts(Ti) = ts(Ti) + tdelta;
        ts(Tc) = ts(Tc) - tdelta;
        if(ts(Ti)<0)
            change = abs(ts(Ti))*0.1;
            ts(Tc) = ts(Tc)-abs(ts(Ti))-change;
            ts(Ti) = change;
        end
        
        if(ts(Tc)<0)
            %ts(Tc) = obj.ts(Tc)*0.1;
            %ts(Ti) = obj.ts(Ti)+obj.ts(Tc)*0.9;
            change = abs(ts(Tc))*0.1;
            ts(Ti) = ts(Ti)-abs(ts(Tc))-change;
            ts(Tc) = change;
            if(ts(Ti)<0)
                ts(Tc) = obj.ts(Tc)*0.1;
                ts(Ti) = obj.ts(Ti)+obj.ts(Tc)*0.9;
            end
        end
        
      %  if(debug), disp(['-- Vsw increasing.  Adjusted dead time Ti=' num2str(Ti) ' by ' num2str(tdelta/sum(ts)*100) '%']); end
    elseif(Vio==0)
        tdelta = (Vmin-Xs(Sir,Xic))/(dxsdt(Sir,Xic));
        % tdelta = min(max(tdelta, -ts(Ti) + delta_DTs), tsmax(Ti) - ts(Ti));
        tdelta = sign(tdelta)*min(abs(tdelta), maxStep);
        if Xic>=3
            if tdelta>0 && Ti<Tc && obj.Converter.Topology.Parser.ONorOFF(Sir,Xic-2)==2
                tdelta = -ts(Ti)/2;
            end
        end
        
        % if tdelta*direction == -1
        %     tdelta = -tdelta;
        % end
        ts(Ti) = ts(Ti) + tdelta;
        ts(Tc) = ts(Tc) - tdelta;
        if(ts(Ti)<0)
            change = abs(ts(Ti))*0.1;
            ts(Tc) = ts(Tc)-abs(ts(Ti))-change;
            ts(Ti) = change;
        end
        if(ts(Tc)<0)
            change = abs(ts(Tc))*0.1;
            ts(Ti) = ts(Ti)-abs(ts(Tc))-change;
            ts(Tc) = change;
            if(ts(Ti)<0)
                ts(Tc) = obj.ts(Tc)*0.1;
                ts(Ti) = obj.ts(Ti)+obj.ts(Tc)*0.9;
            end
        end
        
        
        
       % if(debug), disp(['-- Vsw decreasing.  Adjusted dead time Ti=' num2str(Ti) ' by ' num2str(tdelta/sum(ts)*100) '%']); end
    end
    
    
    %% Check if hard switching is necessary
    % if resonance is carrying voltage in the wrong direction
    %{
            ZVS_dir = (abs(Xs(Sir,Xic-1) - Vmax) < abs(Xs(Sir,Xic-1) - Vmin))*-1 + ...
                (abs(Xs(Sir,Xic-1) - Vmax) > abs(Xs(Sir,Xic-1) - Vmin));
            if(sign(xdot(Sir,1)) ~= sign(ZVS_dir))
                hardSwNecessary = 1;
                % zero deadtime to prevent diode conduction flag
                ts(Tc) = ts(Tc) + ts(Ti);
                ts(Ti) = 0;
            end

            % if we've run out of resonant energy
            if size(xdot,2) > 1
                if(sign(xdot(Sir,end)) ~= sign(xdot(Sir,end-1)))
                    hardSwNecessary = 1;
                end
            end
    %}
    
    % Routine that will lock in Xss if there is a affect in SS
    % waveforms to ensure that there will not be a massive
    % response that will cause the SS solve to go unstable
    
    % If for any state variable the max or min of the previous
    % system is exceeded
    %         if sum(max(Xs,[],2) < max(Xs+dxsdt*tdelta,[],2))>1 || sum(min(Xs,[],2) > min(Xs+dxsdt*tdelta,[],2))>1
    %
    %             % If this previous system is within a tolerance of the
    %             % spacing of the max and min values of the variables
    %             % or a tolerance around the max and min values
    %             % themselves
    %             if sum(abs((max(Xs+dxsdt*tdelta,[],2) - max(Xs,[],2))) > (max(Xs,[],2)-min(Xs,[],2))*0.05)>1 && sum(abs(max(Xs+dxsdt*tdelta,[],2) - max(Xs,[],2)) > abs(max(Xs,[],2)*0.05) )
    %                 keep_SS = true;
    %             end
    %             if sum(abs((min(Xs+dxsdt*tdelta,[],2) - min(Xs,[],2))) > (max(Xs,[],2)-min(Xs,[],2))*0.05)>1 && sum(abs(min(Xs+dxsdt*tdelta,[],2)-min(Xs,[],2)) > abs(min(Xs,[],2)*0.05) )
    %                 keep_SS = true;
    %             end
    %         end
    
    
    tdelta_act = ts(Ti)-obj.ts(Ti);
    
    New_Xs = tdelta_act*dxsdt+Xs;
    obj.Xs = New_Xs;
    
    obj.setts(ts);
    
    %obj.SS_Soln();
    %obj.CorrectXs();
    
    % end % this is the end for the massive overshoot if statement
    return
    
    %end
    
catch ERROR
    %     ERROR;
    rethrow(ERROR);
    %     ts=0; dxsdt=0; hardSwNecessary=0; multcross=0; overresonant=0;
    %     return;
end
end

