function [] = Baxter_adjustDiodeCond(obj)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

% Variable declaration:
Vg = 5;
L = 230e-9;
Cout = 4040e-9;
fs = 2e6;
Ts = 1/fs;
V = 1.8;
Io = 1;
Cp = 3.4874e-10;

Rl =  .01;
dt = Ts/1000;  %50e-9;
Vdr = 5;
Rshunt = 10*270e3;
ron = .05;

u = [Vg Io]';
obj.ts = [Ts*.5-dt dt Ts*.5-dt dt];

ts=obj.ts;

tsmax = ts.*1.2;

As = obj.As;
Bs = obj.Bs;
Cs = obj.Cs;
Ds = obj.Ds;

Xs = obj.Xs;

maxVoErr = 0.05;

debug = 1;

%% Simulate and correct dead time/ duty cycle
%[ Xs] = SS_Soln( As, Bs, ts, u);

% Not other way to find the voltage error of a converter than to be given
% the ouput voltage needed; however, can always be given the postion of
% output cap or something like that.

Voerr = mean(Xs(3,:)) - V;


% The identification of what position these states are in and which ones
% need to be looked at can be automated

% Start by knowing where switching trasnistion will occur:
% Already known by parse.ONorOFF

% From there can determine what states need to be checked and pull these
% measurement nodes from C and D already calcualted

% Remember that FETs volage is from drain to source




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

StateNumbers = obj.Converter.Topology.Parser.StateNumbers;

% Rewrite this to use y instead of x
while (nattempts < 10) && sum(sum(modelError))
	if(debug)
       [xs, t, y] = obj.SS_WF_Reconstruct();
            nattempts
            figure(1)
            ns = size(xs,1);
            for i=1:ns
                subplot(10*ns,1,i*10-9:i*10)
                hold on;
                plot(t,y(StateNumbers(i),:), 'Linewidth', 3); 
                ylabel(obj.getstatenames{i})
                box on

                if(i<ns)
                    set(gca, 'Xticklabel', []);
                else
                    xlabel('t')
                end
            end
            drawnow;


	end

    introduced_Voerr = 0;
    
    if(LSdiode_DT1 || HSdiode_DT1 || hardSw_DT1)
        [tsnew, dxsdtd, hardSwNecessary_DT1] = obj.Baxter_adjustDiodeConduction( Xs, 3, 1, Vg, 0, 1);
        introduced_Voerr = sum(dxsdtd(3,2:end).*(ts-tsnew));
        ts = tsnew;
    end
    
    if(LSdiode_DT2 || HSdiode_DT2 || hardSw_DT2)
        [tsnew, dxsdtd, hardSwNecessary_DT2] = obj.Baxter_adjustDiodeConduction( Xs, 5, 1, Vg, 0, 1);
        introduced_Voerr = introduced_Voerr + sum(dxsdtd(3,2:end).*(ts-tsnew));
        ts = tsnew;
    end
    

    %% compensate Vo for error (original) plus change from dead times
    delta_DTs = max(min(ts)/10, sum(ts)/10000);
    dXs = obj.StateSensitivity('ts', 1, delta_DTs, 3);
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
    obj.ts = ts;
    [ Xs] = obj.SS_Soln();
    
    Voerr = mean(Xs(3,:)) - V;
    LSdiode_DT1 = Xs(1,3) < -10*ron*2;
    HSdiode_DT1 = Xs(1,3) > Vg + 10*ron*2; 

    LSdiode_DT2 = Xs(1,5) < -10*ron*2;  
    HSdiode_DT2 = Xs(1,5) > Vg + 10*ron*2; 
    
    
    [~, ~, ~, mx1] = obj.Baxter_adjustDiodeConduction( Xs, 3, 1, Vg, 0, 1);
    [~, ~, ~, mx2] = obj.Baxter_adjustDiodeConduction( Xs, 5, 1, Vg, 0, 1);
    hardSw_DT1 = (Xs(1,3) > Vg*.01 && ts(2) < tsmax(2) && ~hardSwNecessary_DT1) || mx1;
    hardSw_DT2 = (Xs(1,5) < Vg*.99 && ts(4) < tsmax(4) && ~hardSwNecessary_DT2) || mx2;
    
    modelError = [(abs(Voerr)>maxVoErr) 0 0;
        LSdiode_DT1 HSdiode_DT1 hardSw_DT1;
        LSdiode_DT2 HSdiode_DT2 hardSw_DT2];
    
    nattempts = nattempts + 1;

end

end

