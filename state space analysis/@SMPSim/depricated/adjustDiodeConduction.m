function [ ts, dxsdt, hardSwNecessary, multcross, overresonant] = adjustDiodeConduction(obj, Xs, Xic, Sir, Vmax, Vmin, progBar)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Xic = index to states (Xs) for discrete time interval which (may) contain
% errors due to diode conduction
% Sir = index to states (Xs) for the state in quesiton

    As = obj.As;
    Bs = obj.Bs;
    ts = obj.ts;
    u = obj.u;
    
    tsmax = ts; % JUST FOR TESTING
    
    try
        debug = 0;
        if(debug == 1)
            origts = ts;
        end

        if nargin < 11
            progBar = 1;
        end

        maxStep = sum(ts)/100/(.5+1.5*progBar);

        if ~(Xic>1) || ~(Xic<=size(Xs,2))
            ME = MException('resultisNaN:noSuchVariable', ...
                               'Xi must be between 2 and length(Xs)');
            throw(ME);
        end

        Ti = Xic-1;
        Tc = Ti+1;
        if(Tc == length(ts) + 1)
            Tc = 1;
        end

        dxsdt = zeros(size(Xs));

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

        slopechange = find(diff(sign(xdot(1,:))));

        hardSwNecessary = 0;

        % adjust limits to account for possible I*Ron outside of bounds.
        Vmax = max(Vmax, Xs(Sir, Xic-1));
        Vmin = min(Vmin, Xs(Sir, Xic-1));

        Vmaxcross = find(diff(x(:,Sir)>Vmax));
        Vmincross = find(diff(x(:,Sir)<Vmin));

        multcross = length(Vmaxcross) + length(Vmincross) > 1;
        overresonant = ~isempty(slopechange);
        massiveOvershoot = Xs(Sir,Xic) > Vmax + max(abs([Vmax Vmin 1])) || Xs(Sir,Xic) < Vmin - max(abs([Vmax Vmin 1]));

        if multcross || overresonant || massiveOvershoot
            %% case where deadtime ringing causes multiple diode conduction intervals
            % dead time is way too long, so set it to first crossing (unintelligently) 
            % and try again on next iteration
            tdelta = tsim(max(min([Vmaxcross' Vmincross' slopechange]),2)) - ts(Xic-1);

            %guaranteed to be > 0, so no need to check
            ts(Ti) = ts(Ti) + tdelta;
            ts(Tc) = ts(Tc) - tdelta;
            if(debug), disp('Multicrossing found'); end
            return;
        else

    %     if isempty(Vmaxcross) && isempty(Vmincross) && isempty(slopechange)
    %         % no diode conduction problem

            %% check for ZVS where possible
            delta_DTs = max(min(ts)/100, sum(ts)/10000);
            dXs = obj.StateSensitivity( 'ts', Ti, delta_DTs, Tc);
            dxsdt = (dXs-Xs)/delta_DTs;

            if(xdot(Sir,end) > 0)
                tdelta = (Vmax-Xs(Sir,Xic))/(dxsdt(Sir,Xic));
                tdelta = min(max(tdelta, -ts(Ti) + delta_DTs), tsmax(Ti) - ts(Ti));
                tdelta = sign(tdelta)*min(abs(tdelta), maxStep);
                if(Xs(Sir,Xic)-Vmax < 0 && tdelta < 0)
                   % Slope is positive, and non-ZVS, but partial tells me to reduce dead time
                   % effect will be small (reduced current due to reduced t1)
                   % so tdelta will always saturate and oscillate. Reduce step size 
                   % significantly to prevent oscillating
                   tdelta = sign(tdelta)*min(abs(tdelta), maxStep/10);
                end
                ts(Ti) = ts(Ti) + tdelta;
                ts(Tc) = ts(Tc) - tdelta;
                if(debug), disp(['-- Vsw increasing.  Adjusted dead time Ti=' num2str(Ti) ' by ' num2str(tdelta/sum(ts)*100) '%']); end
            elseif(xdot(Sir,end) < 0)
                tdelta = (Vmin-Xs(Sir,Xic))/(dxsdt(Sir,Xic));
                tdelta = min(max(tdelta, -ts(Ti) + delta_DTs), tsmax(Ti) - ts(Ti));
                tdelta = sign(tdelta)*min(abs(tdelta), maxStep);
                ts(Ti) = ts(Ti) + tdelta;
                ts(Tc) = ts(Tc) - tdelta;
                if(debug), disp(['-- Vsw decreasing.  Adjusted dead time Ti=' num2str(Ti) ' by ' num2str(tdelta/sum(ts)*100) '%']); end
            end


            %% Check if hard switching is necessary
            % if resonance is carrying voltage in the wrong direction
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
            return;    





        end
    catch ERROR
    %     ERROR;
        rethrow(ERROR);
    %     ts=0; dxsdt=0; hardSwNecessary=0; multcross=0; overresonant=0;
    %     return;
    end
end

