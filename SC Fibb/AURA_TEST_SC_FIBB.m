

%     _   _   _  ____    _
%    / \ | | | |/ _  |  / \
%   / _ \| | | | (_| | / _ \
%  / ___ | |_| |> _  |/ ___ \
% /_/   \_\___//_/ |_/_/   \_\



% Current konwn issues use FETs for diodes. (Just have them never turn
% on)

tic
% Place in the filename
filename = 'SC_FIB_AURA.net'; % Place netlist filename here that you want to run

% Can also place component values here to get their voltage and
% current output waveforms such as:

Voltage = {'V1'
    'V2'
    'V3'
    };
Current = {'V1'
    'V2'
    'V3'
    };


%% Determining Coss and ron based on w

%{
w = 0.25;

a =  0.001378;
b = -1;
ron = a*w^b;

p1 = 2.925e-9;
p2 = -2.787e-12;
Coss = p1*w+p2;
%}

%% This is all caluclations to set up the variables need to find the SS
% Solution
Vg = 9.5;


Vb1 = 4;
Vb2 = 4;

Rb = 5e-3;
RL = 5e-3;

ron = 8.5e-3;%8.5e-3;

Coss = 0.9e-9;%.9e-9;

Co = 2e-6; ESRo = 2e-3;
Cfly1 = 5e-6; ESR1 = 3e-3;  % 5V Cap
Cfly2 = 2.5e-6; ESR2 = 3e-3; % 10V Cap
Lc = 100e-9;

fs = 2e6;
Ts = 1/fs;

u = [Vg Vb1 Vb2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]';

modSchemes(:,:,1) = [
    0     1     1     1     0     1     0     1     1     0     0     1     1     1     0     1
    0     1     1     1     0     1     0     1     0     1     0     1     1     1     0     1
    0     1     1     1     0     1     0     1     1     0     0     1     1     1     0     1
    0     1     1     1     0     1     1     0     1     0     0     1     1     1     0     1
    ];

modSchemes(:,:,2) = [
    0     1     1     1     0     1     0     1     0     1     0     1     1     1     0     1
    0     1     1     1     0     1     1     0     0     1     0     1     1     1     0     1
    0     1     1     1     0     1     1     0     1     0     0     1     1     1     0     1
    0     1     1     1     0     1     1     0     0     1     0     1     1     1     0     1
    ];

modSchemes(:,:,3) = [
    0     1     1     1     0     1     1     0     0     1     0     1     1     1     0     1
    0     1     1     0     1     0     1     0     0     1     0     1     1     0     0     1
    0     1     1     1     0     1     1     0     0     1     0     1     1     1     0     1
    1     0     0     1     0     1     1     0     0     1     0     0     1     1     0     1
    ];

modSchemes(:,:,4) = [
    0     1     1     0     1     0     1     0     0     1     0     1     1     0     0     1
    0     1     1     1     0     1     1     0     0     1     0     1     1     0     1     0
    1     0     0     1     0     1     1     0     0     1     0     0     1     1     0     1
    0     1     1     1     0     1     1     0     0     1     1     0     0     1     0     1
    ];

modSchemes(:,:,5) = [
    0     1     1     1     0     1     1     0     0     1     0     1     1     0     1     0
    0     1     1     0     1     0     1     0     0     1     0     1     1     0     1     0
    0     1     1     1     0     1     1     0     0     1     1     0     0     1     0     1
    1     0     0     1     0     1     1     0     0     1     1     0     0     1     0     1
    ];

modSchemes(:,:,6) = [
    0     1     1     0     1     0     1     0     0     1     0     1     1     0     1     0
    0     1     1     1     0     1     1     0     0     1     1     0     0     0     1     0
    1     0     0     1     0     1     1     0     0     1     1     0     0     1     0     1
    0     1     1     1     0     1     1     0     0     1     1     0     0     0     1     0
    ];


fs = 2e6;
Ts = 1/fs;

d = 0.7;

ds = [d, 1-d, d,1-d];

ts = ds/sum(ds)*Ts;


%% There are several variable that must be filled out here:

% You must define u (Input variables in the order of [Independent
% Voltage Soruces MOSFET Forward Votlage (in order of netlist)
% Independent Current Sources]'
%%%% Example u  = [Vg Vfwd1 Vfwd2 Iout]';
u = u;


% Define the inital guess of time intervals. This only defines the
% active switching time. Do not account for diode switching times.
%%%% For example a synchronous Buck converter would be:
%%%%% ts = [Ts*(D-dead) Ts*(dead) Ts*(1-D-dead) Ts*(dead)];
%%%% But a non-synchronous buck covnerter would be
%%%%% ts = [Ts*(D-dead) Ts*(1-(D-dead))];

ts = ts; % The inital guess of time intervals % The inital guess of time intervals



% List all of the numerical components in the netlist file for all
% FETs you must use the syntax used below:
Numerical_Components = {
    'C1' Cfly1
    'C2' Cfly1
    'C3' Cfly2
    'C4' Cfly2
    'C5' Co
    'C6' Co
    'L1' Lc
    'L2' Lc
    'R1' Rb
    'R2' Rb
    'R3' RL
    'R4' RL
    'R5' ESR1
    'R6' ESR2
    'R7' ESR1
    'R8' ESR2
    'R9' ESRo
    'R10' ESRo
    'M1_C' Coss
    'M2_C' Coss
    'M3_C' Coss
    'M4_C' Coss
    'M5_C' Coss
    'M6_C' Coss
    'M7_C' Coss
    'M8_C' Coss
    'M9_C' Coss
    'M10_C' Coss
    'M11_C' Coss
    'M12_C' Coss
    'M13_C' Coss
    'M14_C' Coss
    'M15_C' Coss
    'M16_C' Coss
    };

% List out all char variables in the
Switch_Resistors = {
    'M1_R'
    'M2_R'
    'M3_R'
    'M4_R'
    'M5_R'
    'M6_R'
    'M7_R'
    'M8_R'
    'M9_R'
    'M10_R'
    'M11_R'
    'M12_R'
    'M13_R'
    'M14_R'
    'M15_R'
    'M16_R'
    };
% List of the switch sequency. Organized by: the FETs (column) vs time
% interval (rows) matching Switch_Resistors and ts respectivly

ON = 1;
OFF = 0;
Switch_Sequence = modSchemes(:,:,3);


% List all the resistances of the diodes or FETS when they are on or
% off
SW_OFF = ones(1,16).*10000000;

%SW_ON = [M1_R,M2_R,etc...]
SW_ON = ones(1,16).*ron;
SW = [SW_OFF;SW_ON;SW_ON];

Diode_Forward_Voltage = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]';



%% This creates and runs the parser to set up the converter

parse = NetListParse();
parse.initialize(filename,Voltage,Current);
parse.ABCD();

top = SMPStopology();
top.Parser = parse;

conv = SMPSconverter();
conv.Topology = top;

Order  = 1:length(ts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

conv.ts = ts;
conv.u = u;
conv.order = Order;
conv.Element_Properties = Numerical_Components;
conv.Switch_Resistors = Switch_Resistors;
conv.Switch_Resistor_Values = SW;
conv.Switch_Sequence = Switch_Sequence;
conv.Fwd_Voltage = Diode_Forward_Voltage;


sim = Run_SS_Converter(conv);
toc

return

conv.Element_Properties  = {'C1' C1
    'L1' L1
    'M1_C' M1_C+M1_C
    'M2_C' M2_C
    'R1' R1
    'R2' R2
    'R3' R3
    'M2_C_Resist' M2_C_Resist
    'M1_C_Resist' M1_C_Resist
    };

sim = Run_SS_Converter(conv);

% Example to find the best frequency for the covnerter

X0 = fs_adj;

[X,RESNORM,RESIDUAL,EXITFLAG] = ...
    lsqnonlin(@(x) SS_Error(x, conv),X0);





function [err] = SS_Error(X0, conv)
% This function just accounts for error and runs the function below


Ts = 1/X0;
D = 5/12;
dead = 0.002;
ts = [Ts*(D-dead) Ts*(dead) Ts*(1-D-dead) Ts*(dead)]; % The inital guess of time intervals
conv.ts = ts;
sim = Run_SS_Converter(conv);

[xs, t, y, time_interval] = sim.SS_WF_Reconstruct();

P_in = mean(y(1,:).*y(5,:));

P_out = mean(y(9,:))^2/1;

err = 1-(P_out/P_in);

end



function Plot_Waveforms(sim,statenum,oppstatenum,plotxparam)
%PLOT_WAVEFROMS is a function that plots all of the states
%and the inverse states (V->I or I->V). They will appear in
%figures 100 and 101. (Only input needed is simulation class)




[xs, t, y, time_interval] = sim.SS_WF_Reconstruct();
StateNumbers = sim.Converter.Topology.Parser.StateNumbers;
StateNumbers_Opp = sim.Converter.Topology.Parser.StateNumbers_Opposite;

switch nargin
    case 1
        statenum = 100;
        oppstatenum = 101;
        plotxparam  = 0;
    case 2
        oppstatenum = 101;
        plotxparam = 0;
    case 3
        plotxparam = plotxparam;
end

if plotxparam~=0
    ns = length(plotxparam);
    figure(50)
    for z=1:ns
        ax = subplot(10*ns,1,z*10-9:z*10);
        hold on;
        if plotxparam(z)<0
            plotxparam(z) = abs(plotxparam(z));
            plot(t*10^6,-y(plotxparam(z),:), 'Linewidth', 3);
        else
            plot(t*10^6,y(plotxparam(z),:), 'Linewidth', 3);
        end
        ylabel(plotxparam(z))
        box on
        %ax.YLim = [min(xs(z,:))-abs(0.5*min(xs(z,:))) max(xs(z,:))+abs(0.5*max(xs(z,:)))];
        if(z<ns)
            set(gca, 'Xticklabel', []);
        else
            xlabel('t [\mus]')
        end
    end
end


figure(statenum)
ns = size(xs,1);
for z=1:ns
    ax = subplot(10*ns,1,z*10-9:z*10);
    hold on;
    plot(t,y(StateNumbers(z),:), 'Linewidth', 3);
    ylabel(sim.getstatenames{z})
    box on
    ax.YLim = [min(y(StateNumbers(z),:))-abs(0.5*min(y(StateNumbers(z),:))) max(y(StateNumbers(z),:))+abs(0.5*max(y(StateNumbers(z),:)))];
    if(z<ns)
        set(gca, 'Xticklabel', []);
    else
        xlabel('t(s)')
    end
end
drawnow;


figure(oppstatenum)
ns = size(xs,1);
for z=1:ns
    ax = subplot(10*ns,1,z*10-9:z*10);
    hold on;
    plot(t,y(StateNumbers_Opp(z),:), 'Linewidth', 3);
    ylabel(sim.getstatenames_Opp{z})
    box on
    ax.YLim = [min(y(StateNumbers_Opp(z),:))-abs(0.5*min(y(StateNumbers_Opp(z),:))) max(y(StateNumbers_Opp(z),:))+abs(0.5*max(y(StateNumbers_Opp(z),:)))];
    if(z<ns)
        set(gca, 'Xticklabel', []);
    else
        xlabel('t(s)')
    end
end
drawnow;

end













