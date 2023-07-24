function [ output_args ] = regulate(obj)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

    debug = 1;

    cnstr = obj.converter.topology.constraints;

    As = obj.As;
    Bs = obj.Bs;
    ts = obj.ts;
    u = obj.u;

    %% Simulate and correct dead time/ duty cycle
    [ Xs] = obj.SS_Soln();
    

    for i = 1:size(Xs,2)
        maxcheck(:,i) = cnstr.Cmax*Xs(:,i) + cnstr.Dmax*u - cnstr.Maxlim;
        mincheck(:,i) = cnstr.Cmin*Xs(:,i) + cnstr.Dmin*u + cnstr.Minlim;
        regcheck(:,i) = cnstr.Creg*Xs(:,i) + cnstr.Dreg*u - cnstr.regtarget;
    end
    
    modelError = [
        sum(sum(maxcheck > 0));
        sum(sum(mincheck < 0));
        sum(abs(regcheck) > .1)];
    

    if(sum(modelError(1:2)) ~= 0)
       % Check internal interval waveforms to be sure
       for i = [2 4] %% KNOWN DEADTIME INTERVALS
            [ valid, newt, dist ] = validByInterval(obj, i, Xs);
            diff = newt - obj.ts(i);
%             obj.ts(i) = newt;
%             nextint = max(1,mod(i+1,length(obj.ts)+1));
%             obj.ts(nextint) = obj.ts(nextint) + diff;
            ts = obj.converter.adjustTiming(i,diff);
            [ Xs] = obj.SS_Soln();
            ts(2)

       end
    end

    
%     LSdiode_DT1 = Xs(Sir,3) < -Vfmax;
%     HSdiode_DT1 = Xs(Sir,3) > Vg + Vfmax; 
% 
%     LSdiode_DT2 = Xs(Sir,5) < -max(abs(Xs(3,3:4)))*ronLS*2;  
%     HSdiode_DT2 = Xs(Sir,5) > Vg + max(abs(Xs(3,3:4)))*ronHS*2; 
% 
%     hardSwNecessary_DT1 = 0;
%     hardSwNecessary_DT2 = 0;
%     hardSw_DT1 = Xs(1,3) > Vg*maxVswErr && ts(2) < tsmax(2) && ~hardSwNecessary_DT1;
%     hardSw_DT2 = Xs(1,5) < Vg*(1-maxVswErr) && ts(4) < tsmax(4) && ~hardSwNecessary_DT2;
% 
%     modelError = [(abs(Voerr)>maxVoErr) 0 0;
%         LSdiode_DT1 HSdiode_DT1 hardSw_DT1;
%         LSdiode_DT2 HSdiode_DT2 hardSw_DT2];

    nattempts = 0;
    maxNattempts = 15;

    if(debug)
        figure(2);
        hold off;
        set(gcf,'position',[1921,0,1919,1200])
    end

    while (nattempts < maxNattempts) && sum(sum(modelError))

        if(debug)
            figure(obj.debugFigNo);
            set(gcf,'position',[1921,0,1919,1200])
            
            obj.plotAllStates(obj.debugFigNo, 0);
            disp('------------------------------------------------------');
            disp(['Iteration number ' num2str(nattempts) ' of ' num2str(maxNattempts)]); 
        end

        introduced_regerr = 0;
%         hardSwNecessary_DT2 = 0;
%         hardSwNecessary_DT1 = 0;

        %% calculate sensitivity based on original ts
        delta_DTs = max(min(ts)/10, sum(ts)/10000);
        dXs = obj.StateSensitivity( 'ts', 1, delta_DTs, 3);
        dxsdt = (dXs - Xs)/delta_DTs;

        mx1 = 0;
        mx2 = 0;
        tsnew1 = ts;
        tsnew2 = ts;
        
        dt = -ts(2)*99/100;
%         [dt] = obj.converter.adjustTiming(2, -ts(2)*99/100);
        %% Adjust dead times
        if(LSdiode_DT1 || HSdiode_DT1 || hardSw_DT1)
            [tsnew1, dxsdtd1, hardSwNecessary_DT1, mx1] = adjustDiodeConduction( Xs, As, Bs, u, ts, tsmax, 3, 1, Vg, 0, nattempts/maxNattempts);
    %         introduced_Voerr = sum(dxsdtd1(3,2:end).*(ts-tsnew));
            introduced_Voerr = mean(dxsdtd1(3,2:end))*(tsnew1(2)-ts(2));
        end

        if(LSdiode_DT2 || HSdiode_DT2 || hardSw_DT2)
            [tsnew2, dxsdtd2, hardSwNecessary_DT2, mx2] = adjustDiodeConduction( Xs, As, Bs, u, ts, tsmax, 5, 1, Vg, 0, nattempts/maxNattempts);
    %         introduced_Voerr = introduced_Voerr + sum(dxsdtd2(3,2:end).*(ts-tsnew));
            delta_t4 = tsnew2(4)-ts(4);
            introduced_Voerr = introduced_Voerr + mean(dxsdtd2(3,2:end))*(tsnew2(4)-ts(4));
        end

        ts = [tsnew2(1) tsnew1(2) tsnew1(3) tsnew2(4)];


        %% compensate Vo for error (original) plus change from dead times
        dt = -(Voerr+introduced_Voerr)/mean(dxsdt(3,2:end));

        if(mx1 || mx2 || Xs(1,3) < -3 || Xs(1,5) < -3 || Xs(1,3) > Vg*2 || Xs(1,5) > Vg*2)
            % dead times way too long (multiple resonant periods).  So, don't
            % mess with duty cycle until they're sorted.
            dt = 0;
        end
        %% Check if full soft switching is possible (by IL) but cannot simultaneously regulate

        if (LSdiode_DT2 || HSdiode_DT2 || hardSw_DT2) && ~(LSdiode_DT1 || HSdiode_DT1 || hardSw_DT1)
            % meant to increase V5 by
            delta_v5_dt = delta_t4*dxsdtd2(1,5);
            % to keep output voltage, it will further change by
            delta_v5_D = dxsdt(1,5)*dt;

            if( abs(delta_v5_D) > abs(delta_v5_dt) && sign(delta_v5_D) + sign(delta_v5_dt) == 0)
               %can't get ZVS without losing regulation -- going to have to hard switch
               new_dt4 =  dxsdt(1,5)*Voerr/(-dxsdtd2(1,5)*mean(dxsdt(3,2:end)) - dxsdt(1,5)*mean(-dxsdtd2(3,2:end)));
               new_dt4 = max(0,new_dt4);
               xchange = new_dt4 + delta_t4;
               ts(4) = ts(4) + xchange;
               ts(1) = ts(1) - xchange;
               hardSwNecessary_DT2 = 1;
               introduced_Voerr =  mean(dxsdtd2(3,2:end))*new_dt4;
               dt = -(Voerr+introduced_Voerr)/mean(dxsdt(3,2:end));
               if(debug), disp('-----Hard switching seems to be necessary.  Trying to minimize non-ZVS voltage'); end
            end
        end

    %     dt = sign(dt)*min(abs(dt), min(.95*ts(1), .95*ts(3)));

        %% Adjust D*Ts 
        if(ts(3) - dt <0)
            dt = ts(3);
        elseif(ts(1) + dt < 0)
            dt = -ts(1);
        end

        if(debug)
            disp(['Voerr = ' num2str(Voerr) ' and introduced_Voerr = ' num2str(introduced_Voerr)]);
            if(Voerr+introduced_Voerr > 0)
                disp(['- Vo too large.  Adjusted DTs by ' num2str(dt/sum(ts)*100) '%']); 
            elseif(Voerr+introduced_Voerr < 0)
                disp(['- Vo too small.  Adjusted DTs by ' num2str(dt/sum(ts)*100) '%']); 
            end
        end


        ts(1) = ts(1) + dt;
        ts(3) = ts(3) - dt;

        if(ts(3) > tsmax(3))
            ts(1) = ts(1) - tsmax(3) + ts(3);
            ts(3) = tsmax(3);
            if(debug), disp(['clipped due to t3 exceeding 90%']); end
        elseif(ts(1) > tsmax(1))
            ts(3) = ts(3) - tsmax(1) + ts(1);
            ts(1) = tsmax(1);
            if(debug), disp(['----- clipped due to t1 exceeding 90%']); end
        end


        %% Recompute and reevaluate
        [ Xs] = obj.SS_Soln();

%         Voerr = mean(Xs(3,2:end)) - V;
%         LSdiode_DT1 = Xs(1,3) < -max(abs(Xs(2,2:3)))*ronLS*2;
%         HSdiode_DT1 = Xs(1,3) > Vg + max(abs(Xs(2,2:3)))*ronHS*2; 
% 
%         LSdiode_DT2 = Xs(1,5) < -max(abs(Xs(2,3:4)))*ronLS*2;  
%         HSdiode_DT2 = Xs(1,5) > Vg + max(abs(Xs(2,3:4)))*ronHS*2; 
    
        for i = 1:size(Xs,2)
            maxcheck(:,i) = cnstr.Cmax*Xs(:,i) + cnstr.Dmax*u - cnstr.Maxlim;
            mincheck(:,i) = cnstr.Cmin*Xs(:,i) + cnstr.Dmin*u + cnstr.Minlim;
            regcheck(:,i) = cnstr.Creg*Xs(:,i) + cnstr.Dreg*u - cnstr.regtarget;
        end

        modelError = [
            sum(sum(maxcheck > 0));
            sum(sum(maxcheck < 0));
            sum(abs(regcheck) > .1)];


%         %Doesn't change the times, but double-check for multicrossing before
%         %allowing a solution
%         [~, ~, ~, mx1] = adjustDiodeConduction( Xs, As, Bs, u, ts, tsmax, 3, 1, Vg, 0, 1);
%         [~, ~, ~, mx2] = adjustDiodeConduction( Xs, As, Bs, u, ts, tsmax, 5, 1, Vg, 0, 1);


%         hardSw_DT1 = (Xs(1,3) > Vg*maxVswErr && ts(2) < tsmax(2) && ~hardSwNecessary_DT1) || mx1;
%         hardSw_DT2 = (Xs(1,5) < Vg*(1-maxVswErr) && ts(4) < tsmax(4) && ~hardSwNecessary_DT2) || mx2;
% 
%         modelError = [(abs(Voerr)>maxVoErr) 0 0;
%             LSdiode_DT1 HSdiode_DT1 hardSw_DT1;
%             LSdiode_DT2 HSdiode_DT2 hardSw_DT2];

        nattempts = nattempts + 1;

        if(nattempts > maxNattempts-5) 
            %likely going to fail; de-prioritize ZVs to get regulation
%             hardSw_DT1 =  mx1;
%             hardSw_DT2 =  mx2;
            if(debug), disp(['----- Not looking good!; give up on ZVS']); end
        end





    end

end

