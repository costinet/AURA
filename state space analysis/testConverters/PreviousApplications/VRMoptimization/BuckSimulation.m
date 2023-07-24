function [ eta, Xs, ts, xst, ys, t, exitflag ] = BuckSimulation( inductor, HSFET, LSFET, fs, Vg, V, Io, Cout, optOpts )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if(nargin == 0)
    BuckOptimization;
end

if(fs<0)
    eta = -1;
    Xs = 0;
    ts = 0;
    xst = 0;
    ys = 0;
    t = 0;
    exitflag = -1;
    return;
end

maxVoErr = .01;
maxVswErr = .01;

debug = 1;

% Vg = 5;
L = inductor.L;
% Cout = 4040e-9;
% fs = 2e6;
Ts = 1/fs;
% V = 1.8;
% Io = 1;
CpHS = HSFET.Cp;
CpLS = LSFET.Cp;

Rl =  inductor.Rl;
dt = Ts/4;%50e-9;
Rshunt = 10*270e3;
ronHS = HSFET.ron;
ronLS = LSFET.ron;


dt = min(pi*sqrt(inductor.L*(HSFET.Cp+LSFET.Cp))/2, Ts/10);
ts = [V/Vg*Ts-dt dt (1-V/Vg)*Ts-dt dt];
% tsmax = [.9 .25 .9 .25]*Ts;
tsmax = [.9*Ts 3.9791e-08 .9*Ts 3.9791e-08];

%% Waveform reconstruction - Nonsingular Case
% x = [Vp Il Vo]
u = [Vg Io]';


A1 = [-1/ronHS -1 0; 1 -Rl -1; 0 1 0];
A2 = [0 -1 0; 1 -Rl -1; 0 1 -1/Rshunt];
A3 = [-1/ronLS -1 0; 1 -Rl -1; 0 1 0];
A4 = [-1/Rshunt 0 0; 1 -Rl -1; 0 1 -1/Rshunt];

B1 = [1/ronHS 0; 0 0; 0 -1];
B2 = [0 0; 0 0; 0 -1];
B3 = [0 0; 0 0; 0 -1];
B4 = [0 0; 0 0; 0 -1];

K = [CpHS+CpLS 0 0 ; 0 L 0 ; 0 0 Cout];

    

As = cat(3, A1, A2, A3, A2);
Bs = cat(3,B1, B2, B3, B2);

     
for i = 1:size(As,3)
    As(:,:,i) = K^-1*As(:,:,i);
    Bs(:,:,i) = K^-1*Bs(:,:,i);
end


%% Simulate and correct dead time/ duty cycle
[ Xs] = SS_Soln( As, Bs, ts, u);
Voerr = mean(Xs(3,2:end)) - V;
LSdiode_DT1 = Xs(1,3) < -max(abs(Xs(2,2:3)))*ronLS*2;
HSdiode_DT1 = Xs(1,3) > Vg + max(abs(Xs(2,2:3)))*ronHS*2; 

LSdiode_DT2 = Xs(1,5) < -max(abs(Xs(2,3:4)))*ronLS*2;  
HSdiode_DT2 = Xs(1,5) > Vg + max(abs(Xs(2,3:4)))*ronHS*2; 

hardSwNecessary_DT1 = 0;
hardSwNecessary_DT2 = 0;
hardSw_DT1 = Xs(1,3) > Vg*maxVswErr && ts(2) < tsmax(2) && ~hardSwNecessary_DT1;
hardSw_DT2 = Xs(1,5) < Vg*(1-maxVswErr) && ts(4) < tsmax(4) && ~hardSwNecessary_DT2;

modelError = [(abs(Voerr)>maxVoErr) 0 0;
    LSdiode_DT1 HSdiode_DT1 hardSw_DT1;
    LSdiode_DT2 HSdiode_DT2 hardSw_DT2];
 
nattempts = 0;
maxNattempts = 15;

if(debug)
    figure(2);
    hold off;
    set(gcf,'position',[1921,0,1919,1200])
end

while (nattempts < maxNattempts) && sum(sum(modelError))
    
	if(debug)
        [ ys, t ] = SS_WF_Reconstruct( Xs, As, Bs, ts, u );
        disp(modelError);
        disp(ts/Ts);
        lines = findobj(gca,'Type','Line');
        for i = 1:numel(lines)
          lines(i).LineWidth = 1;
          lines(i).LineStyle = ':';
        end
        figure(2)
        plot(t,ys, 'linewidth',2);
        hold on;
        plot(t, V*ones(size(t)), '-.k', 'linewidth',3);
        ylims = [min(min(ys)) max(max(ys))];
        ylim(ylims)
        for i = 1:length(ts)
            plot(sum(ts(1:i))*ones(1,2), ylims, ':r');
        end
        disp('------------------------------------------------------');
        disp(['Iteration number ' num2str(nattempts) ' of ' num2str(maxNattempts)]); 
        pause
	end

    introduced_Voerr = 0;
    hardSwNecessary_DT2 = 0;
    hardSwNecessary_DT1 = 0;
    
    %% calculate sensitivity based on original ts
    delta_DTs = max(min(ts)/10, sum(ts)/10000);
    dXs = StateSensitivity( As, Bs, ts, u, 'ts', 1, delta_DTs, 3);
    dxsdt = (dXs - Xs)/delta_DTs;
    
    mx1 = 0;
    mx2 = 0;
    tsnew1 = ts;
    tsnew2 = ts;
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
    [ Xs] = SS_Soln( As, Bs, ts, u);
    
    Voerr = mean(Xs(3,2:end)) - V;
    LSdiode_DT1 = Xs(1,3) < -max(abs(Xs(2,2:3)))*ronLS*2;
    HSdiode_DT1 = Xs(1,3) > Vg + max(abs(Xs(2,2:3)))*ronHS*2; 

    LSdiode_DT2 = Xs(1,5) < -max(abs(Xs(2,3:4)))*ronLS*2;  
    HSdiode_DT2 = Xs(1,5) > Vg + max(abs(Xs(2,3:4)))*ronHS*2; 

    
    %Doesn't change the times, but double-check for multicrossing before
    %allowing a solution
    [~, ~, ~, mx1] = adjustDiodeConduction( Xs, As, Bs, u, ts, tsmax, 3, 1, Vg, 0, 1);
    [~, ~, ~, mx2] = adjustDiodeConduction( Xs, As, Bs, u, ts, tsmax, 5, 1, Vg, 0, 1);
    
    
    hardSw_DT1 = (Xs(1,3) > Vg*maxVswErr && ts(2) < tsmax(2) && ~hardSwNecessary_DT1) || mx1;
    hardSw_DT2 = (Xs(1,5) < Vg*(1-maxVswErr) && ts(4) < tsmax(4) && ~hardSwNecessary_DT2) || mx2;
    
    modelError = [(abs(Voerr)>maxVoErr) 0 0;
        LSdiode_DT1 HSdiode_DT1 hardSw_DT1;
        LSdiode_DT2 HSdiode_DT2 hardSw_DT2];
    
	nattempts = nattempts + 1;
    
	if(nattempts > maxNattempts-5) 
        %likely going to fail; de-prioritize ZVs to get regulation
        hardSw_DT1 =  mx1;
        hardSw_DT2 =  mx2;
        if(debug), disp(['----- Not looking good!; give up on ZVS']); end
	end
    

    


end




Cp = (CpHS + CpLS);
uCr = CpHS/(CpHS + CpLS);
lCr = CpLS/(CpHS + CpLS);

% x = [Vp Il Vo]
% y = [vo ig];
u = [Vg Io]';


% C1 = [0 0 1; -1/ronHS/Cp*lCr, -1/Cp*lCr + 1, 0];
% D1 = [0 0; 1/ronHS/Cp*lCr, 0];
% 
% C2 = [0 0 1; 0, -1/Cp*uCr, 0];
% D2 = [0 0; 0, 0];
% 
% C3 = [0 0 1; -1/ronLS/Cp*uCr, -1/Cp*uCr, 0];
% D3 = [0 0; 0, 0];

C1 = [0 0 1; (As(1,:,1)*lCr*Cp + [0 1 0])];
D1 = [0 0; Bs(1,:,1)*lCr*Cp];

C2 = [0 0 1; -As(1,:,2)*uCr*Cp];
D2 = [0 0; -Bs(1,:,2)*uCr*Cp];

C3 = [0 0 1; -As(1,:,3)*uCr*Cp];
D3 = [0 0; -Bs(1,:,3)*uCr*Cp];

Cs = cat(3, C1, C2, C3, C2);
Ds = cat(3, D1, D2, D3, D2);

[ xst, t, ys ] = SS_WF_Reconstruct( Xs, As, Bs, ts, u, Cs, Ds );

eta = mean(Io*ys(1,:))/mean(Vg*ys(2,:));

if( eta>1 || eta<0 || Io<=0 || mean(ys(2,:)) <=0)
    eta = -1;
end

if(debug)
    disp(modelError);
    disp(ts);
    lines = findobj(gca,'Type','Line');
    for i = 1:numel(lines)
      lines(i).LineWidth = 1;
      lines(i).LineStyle = ':';
    end
    figure(2)
    plot(t,xst, 'linewidth',8);
    hold on;
    ylims = [min(min(xst)) max(max(xst))];
    ylim(ylims)
    for i = 1:length(ts)
        plot(sum(ts(1:i))*ones(1,2), ylims, ':r');
    end
    pause
end


exitflag.values = [modelError; 0, sum(ts>=0)~=length(ts) nattempts];
exitflag.names = {'Vo Error', 'N/A', 'N/A'; 
        'LS Diode Conduction in DT1', 'HS Diode Conduction in DT1', 'Unneccessary Hard Switching in DT1';
        'LS Diode Conduction in DT2', 'HS Diode Conduction in DT2', 'Unneccessary Hard Switching in DT2';
        'N/A', 'Negative timing intervals', 'number of iterations'};    
    

% plot(t,ys)
% hold on;
% ylims = ylim;
% for i = 1:length(ts)
%     plot(sum(ts(1:i))*ones(1,2), ylims, ':r');
% end
% legend('V_p', 'I_L', 'V_o');








end

