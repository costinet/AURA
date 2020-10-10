

%     _   _   _  ____    _    
%    / \ | | | |/ _  |  / \   
%   / _ \| | | | (_| | / _ \  
%  / ___ | |_| |> _  |/ ___ \ 
% /_/   \_\___//_/ |_/_/   \_\



% Current konwn issues use FETs for diodes. (Just have them never turn
% on)

tic
% Place in the filename
filename = 'HDSC_AURA.net'; % Place netlist filename here that you want to run

% Can also place component values here to get their voltage and
% current output waveforms such as: 

Voltage = {'V1'
    };
Current = {'V1'
    };


%% This is all caluclations to set up the variables need to find the SS
% Solution

Vg = 24;
V = 5;
Io = 10;
u = [Vg 1 1 1 1 1 1 1 1]';
L1 = 2*72e-9;
RL = (0.007-2.65e-3 + 2*.35e-3);
Da = 1;
M = 5/6;

fs = 500e3;
Ts = 1/fs;

dt = 1/fs/Ts;

C4 = 100e-6; % Cout
C1 = 13.38e-6;
% C1 = 1e-6;
C2 = 12.49e-6;
C3 = 10.27e-6;

ESRx = 3e-3;
R3 = 1.1e-3+ESRx;
R2 = 2.3e-3+ESRx;
R4 = 2.3e-3+ESRx;
ron = 2.5e-3;
R1 = 0.5;
R5 = 0.001;
Lp = 0e-9;
Lp1 = Lp;
Lp2 = 2*Lp;
Lp3 = 3*Lp;


kc = 1;
Coss = kc*2.5e-9;
Coss6 = kc*2.26e-9;
Coss12 = kc*1.99e-9;
Qg = 18e-9;

Phase1a = [1 0 1 1 0 0 0 1];
Phase1b = [1 0 0 1 0 0 0 1];
Phase2a = 1 - Phase1a;
Phase2b = Phase2a; Phase2b(6) = 0;
Phase3 = [0 0 0 1 1 0 1 1];
swseq = [Phase1a; Phase1b; Phase3;  Phase2a; Phase2b;  Phase3];


M1_C = Coss6; % CHS
M2_C = Coss12; % LHS
M3_C = Coss12; % CHS
M4_C = Coss6; % LHS
M5_C = Coss6; % CHS
M6_C = Coss6; % LHS
M7_C = Coss6; % CHS
M8_C = Coss6; % LHS

Order = [1 2 3 4 5 6]; % The order that the states must go in after being parsed
Order = [1 2 3 4];
SW_OFF = ones(1,8).*10000000;

SW_ON = ones(1,8).*2.5e-3;

SW = [SW_OFF;SW_ON];

M1 = [1
    1
    1
    0
    0
    1];

M4 = M1;

M2 = [0
    0
    1
    1
    1
    1];

M3 = M2;

M5 = [0
    0
    0
    1
    0
    0];

M6 = [1
    1
    0
    0
    0
    0];

M7 = [0
    0
    0
    1
    1
    0];

M8 = [1
    0
    0
    0
    0
    0];

Binary_for_DAB = [M1 M2 M3 M4 M5 M6 M7 M8]; % But its really for HDSC

Binary_for_DAB(5,:) = [];
Binary_for_DAB(2,:) = [];

ts = [Da*Ts/2*M, (1-Da)*Ts/2*M, Ts/2*(1-M), Da*Ts/2*M, (1-Da)*Ts/2*M, Ts/2*(1-M)];
ts(5) = [];
ts(2) = [];






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
            'C1' C1
            'C2' C2
            'C3' C3
            'C4' C4
            'L1' L1
            'R1' R1
            'R2' R2
            'R3' R3
            'R4' R4
            'R5' R5
            'M1_C' M1_C
            'M2_C' M2_C
            'M3_C' M3_C
            'M4_C' M4_C
            'M5_C' M5_C
            'M6_C' M6_C
            'M7_C' M7_C
            'M8_C' M8_C
            };

% List out all char variables in the 
        Switch_Resistors = {'M1_R'
            'M2_R'
            'M3_R'
            'M4_R'
            'M5_R'
            'M6_R'
            'M7_R'
            'M8_R'};
% List of the switch sequency. Organized by: the FETs (column) vs time
% interval (rows) matching Switch_Resistors and ts respectivly

ON = 1;
OFF = 0;
Switch_Sequence = Binary_for_DAB;


% List all the resistances of the diodes or FETS when they are on or
% off
SW_OFF = ones(1,8).*10000000;

%SW_ON = [M1_R,M2_R];
SW_ON = ones(1,8).*2.5e-3;
SW = [SW_OFF;SW_ON;SW_ON];

Diode_Forward_Voltage = [1 1 1 1 1 1 1 1]';



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













