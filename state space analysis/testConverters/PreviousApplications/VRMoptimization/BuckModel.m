clear all; clc;
% I_L V_DS V_GS V_OUT

useCw = 0;

debug = 0;

Vg = 5;
L = 230e-9;
Cout = 4040e-9;
fs = 2e6;
Ts = 1/fs;
V = 1.8;
Io = 1;
Cp = 3.4874e-10;

Rl =  .01;
dt = Ts/4;%50e-9;
Vdr = 5;
Rshunt = 10*270e3;
ron = .05;

maxVoErr = .01;

ts = [Ts*.5-dt dt Ts*.5-dt dt];
tsmax = ts;

%% Waveform reconstruction - Nonsingular Case
% x = [Vp Il Vo]
u = [Vg Io]';


A1 = [-1/ron -1 0; 1 -Rl -1; 0 1 0];
A2 = [0 -1 0; 1 -Rl -1; 0 1 -1/Rshunt];
A3 = [-1/ron -1 0; 1 -Rl -1; 0 1 0];
A4 = [-1/Rshunt 0 0; 1 -Rl -1; 0 1 -1/Rshunt];

B1 = [1/ron 0; 0 0; 0 -1];
B2 = [0 0; 0 0; 0 -1];
B3 = [0 0; 0 0; 0 -1];
B4 = [0 0; 0 0; 0 -1];

K = [Cp 0 0 ; 0 L 0 ; 0 0 Cout];

    

As = cat(3, A1, A2, A3, A2);
Bs = cat(3,B1, B2, B3, B2);

     
for i = 1:size(As,3)
    As(:,:,i) = K^-1*As(:,:,i);
    Bs(:,:,i) = K^-1*Bs(:,:,i);
end


%% Simulate and correct dead time/ duty cycle
[ Xs] = SS_Soln( As, Bs, ts, u);
Voerr = mean(Xs(3,:)) - V;
LSdiode_DT1 = Xs(1,3) < -10*ron*2;
HSdiode_DT1 = Xs(1,3) > Vg + 10*ron*2; 

LSdiode_DT2 = Xs(1,5) < -10*ron*2;  
HSdiode_DT2 = Xs(1,5) > Vg + 10*ron*2; 

hardSwNecessary_DT1 = 0;
hardSwNecessary_DT2 = 0;
hardSw_DT1 = Xs(1,3) > Vg*.01 && ts(2) < tsmax(2) && ~hardSwNecessary_DT1;
hardSw_DT2 = Xs(1,5) < Vg*.99 && ts(4) < tsmax(4) && ~hardSwNecessary_DT2;

modelError = [(abs(Voerr)>maxVoErr) 0 0;
    LSdiode_DT1 HSdiode_DT1 hardSw_DT1;
    LSdiode_DT2 HSdiode_DT2 hardSw_DT2];
 
nattempts = 0;

while (nattempts < 100) && sum(sum(modelError))
    
	if(debug)
        [ ys, t ] = SS_WF_Reconstruct( Xs, As, Bs, ts, u );
        disp(modelError);
        disp(ts);
        plot(t,ys);
        hold on;
        ylims = ylim;
        for i = 1:length(ts)
            plot(sum(ts(1:i))*ones(1,2), ylims, ':r');
        end
%         pause
	end

    introduced_Voerr = 0;
    
    if(LSdiode_DT1 || HSdiode_DT1 || hardSw_DT1)
        [tsnew, dxsdtd, hardSwNecessary_DT1] = adjustDiodeConduction( Xs, As, Bs, u, ts, tsmax, 3, 1, Vg, 0);
        introduced_Voerr = sum(dxsdtd(3,2:end).*(ts-tsnew));
        ts = tsnew;
    end
    
    if(LSdiode_DT2 || HSdiode_DT2 || hardSw_DT2)
        [tsnew, dxsdtd, hardSwNecessary_DT2] = adjustDiodeConduction( Xs, As, Bs, u, ts, tsmax, 5, 1, Vg, 0);
        introduced_Voerr = introduced_Voerr + sum(dxsdtd(3,2:end).*(ts-tsnew));
        ts = tsnew;
    end
    

    %% compensate Vo for error (original) plus change from dead times
    delta_DTs = max(min(ts)/10, sum(ts)/10000);
    dXs = StateSensitivity( As, Bs, ts, u, 'ts', 1, delta_DTs, 3);
    dxsdt = (Xs - dXs)/delta_DTs;
    dt = (Voerr+introduced_Voerr)/mean(dxsdt(3,:));

    if(ts(3) - dt <0)
        dt = ts(3);
    elseif(ts(1) + dt < 0)
        dt = -ts(1);
    end
    dt = dt*.5;

    ts(1) = ts(1) + dt;
    ts(3) = ts(3) - dt;
    
    %% Recompute and reevaluate
    [ Xs] = SS_Soln( As, Bs, ts, u);
    
    Voerr = mean(Xs(3,:)) - V;
    LSdiode_DT1 = Xs(1,3) < -10*ron*2;
    HSdiode_DT1 = Xs(1,3) > Vg + 10*ron*2; 

    LSdiode_DT2 = Xs(1,5) < -10*ron*2;  
    HSdiode_DT2 = Xs(1,5) > Vg + 10*ron*2; 
    
    
    [~, ~, ~, mx1] = adjustDiodeConduction( Xs, As, Bs, u, ts, tsmax, 3, 1, Vg, 0);
    [~, ~, ~, mx2] = adjustDiodeConduction( Xs, As, Bs, u, ts, tsmax, 5, 1, Vg, 0);
    hardSw_DT1 = (Xs(1,3) > Vg*.01 && ts(2) < tsmax(2) && ~hardSwNecessary_DT1) || mx1;
    hardSw_DT2 = (Xs(1,5) < Vg*.99 && ts(4) < tsmax(4) && ~hardSwNecessary_DT2) || mx2;
    
    modelError = [(abs(Voerr)>maxVoErr) 0 0;
        LSdiode_DT1 HSdiode_DT1 hardSw_DT1;
        LSdiode_DT2 HSdiode_DT2 hardSw_DT2];
    
    nattempts = nattempts + 1;

end

[ ys, t ] = SS_WF_Reconstruct( Xs, As, Bs, ts, u );

plot(t,ys)
hold on;
ylims = ylim;
for i = 1:length(ts)
    plot(sum(ts(1:i))*ones(1,2), ylims, ':r');
end
legend('V_p', 'I_L', 'V_o');

% %% Waveform reconstruction - Singular Case
% % x = [Vp Il Vo]
% 
% A2 = [0 -1 0; 1 -Rl -1; 0 1 0];
% A4 = [0 0 0; 1 -Rl -1; 0 1 0];
% 
% As = cat(3, A1, A2, A3, A4);
%      
% for i = 1:size(As,3)
%     As(:,:,i) = K^-1*As(:,:,i);
% end
% 
% [ Xs2] = SS_Soln( As, Bs, ts, u)
% [ ys, t ] = SS_WF_Reconstruct( Xs2, As, Bs, ts, u );
% 
% plot(t,ys)
% Xs-Xs2





