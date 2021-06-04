

%     _   _   _  ____    _    
%    / \ | | | |/ _  |  / \   
%   / _ \| | | | (_| | / _ \  
%  / ___ | |_| |> _  |/ ___ \ 
% /_/   \_\___//_/ |_/_/   \_\



% Current konwn issues use FETs for diodes. (Just have them never turn
% on)

tic
% Place in the filename
filename = 'Buck_Qual.net'; % Place netlist filename here that you want to run

% Can also place component values here to get their voltage and
% current output waveforms such as: 

Voltage = {'V1'
    };
Current = {'V1'
    };


%% This is all caluclations to set up the variables need to find the SS
% Solution

D = 5/12;
fs_adj = 200e3;
Ts = 1/fs_adj;
dead = 0.002;
Vg = 12;
fs = 0.2e6;
Ts = 1/fs;
dt = 0.02/fs;
L1 = 6e-6;
%L1 = 0.5e-6;
r_on = 0.0015;
Cds = 1620e-12;
C1 = (4*100+680)*10^-6; % Output Cap
R1 = 1;
R2 = 0.16;
R3 = 0.05;
[M1_R,M2_R] = deal(r_on);
[M1_C,M2_C] = deal(Cds);
M2_C_Resist = 10000;
M1_C_Resist = 10000;


%% There are several variable that must be filled out here:

% You must define u (Input variables in the order of [Independent
% Voltage Soruces MOSFET Forward Votlage (in order of netlist)
% Independent Current Sources]'
%%%% Example u  = [Vg Vfwd1 Vfwd2 Iout]';
u = [Vg 1 1]';


% Define the inital guess of time intervals. This only defines the
% active switching time. Do not account for diode switching times.
%%%% For example a synchronous Buck converter would be: 
%%%%% ts = [Ts*(D-dead) Ts*(dead) Ts*(1-D-dead) Ts*(dead)];
%%%% But a non-synchronous buck covnerter would be
%%%%% ts = [Ts*(D-dead) Ts*(1-(D-dead))];

ts = [Ts*(D-dead) Ts*(dead) Ts*(1-D-dead) Ts*(dead)]; % The inital guess of time intervals



% List all of the numerical components in the netlist file for all
% FETs you must use the syntax used below:
Numerical_Components = {'C1' C1
    'L1' L1
    'M1_C' M1_C
    'M2_C' M2_C
    'R1' R1
    'R2' R2
    'R3' R3
    'M1_R' M1_R
    'M2_R' M2_R
    };

% List out all char variables in the 
Switch_Resistors = {'M1_R'
    'M2_R'};

% List of the switch sequency. Organized by: the FETs (column) vs time
% interval (rows) matching Switch_Resistors and ts respectivly

ON = 1;
OFF = 0;
Switch_Sequence = [
    ON OFF
    OFF OFF
    OFF ON
    OFF OFF
    ];


% List all the resistances of the diodes or FETS when they are on or
% off
SW_OFF = ones(1,2).*10000000;

SW_ON = [M1_R,M2_R];

SW = [SW_OFF;SW_ON;SW_ON];


Diode_Forward_Voltage = [1 1]';



%% This creates and runs the parser to set up the converter

parse = NetListParse(); 
parse.initialize(filename,Voltage,Current);
parse.Component_Values = Numerical_Components; %%% Added to test numerical stuff in code
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













