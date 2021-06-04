

%     _   _   _  ____    _    
%    / \ | | | |/ _  |  / \   
%   / _ \| | | | (_| | / _ \  
%  / ___ | |_| |> _  |/ ___ \ 
% /_/   \_\___//_/ |_/_/   \_\



% Current konwn issues use FETs for diodes. (Just have them never turn
% on)

tic
% Place in the filename
filename = 'DAB_Resistors_Cap.net'; % Place netlist filename here that you want to run

% Can also place component values here to get their voltage and
% current output waveforms such as: 

Voltage = {'V1'
    };
Current = {'V1'
    };


%% This is all caluclations to set up the variables need to find the SS
% Solution

 Vg = 48;
        L3 = 1e-6;
        L2 = 1296e-6;
        L1 = 76e-6;
        C1 = (4*100+680)*10^-6; %Cout
        C2 = 100*10^-6;
        fs = 0.2e6;
        Ts = 1/fs;
        V = 1.2;
        Io = 10; % was 1
        M1_C = 110e-12; % CHS
        M2_C = 110e-12; % LHS
        M3_C = 110e-12; % CHS
        M4_C = 110e-12; % LHS
        M5_C = 1620e-12; % CHS
        M6_C = 1620e-12; % LHS
        M7_C = 1620e-12; % CHS
        M8_C = 1620e-12; % LHS
        D1_C = 1620e-12; % LHS
        
        R3 = 0.001;
        
        R2 = 4.57535; % the inductor resistance % 1 ohm on the primary 1Ohm on the secondary
        dt = Ts/1000;%5e-10;
        
        Order = [1 2 3 4 5 6 7 8]; % The order that the states must go in after being parsed
        
        chooses = 9;
        switch chooses
            
            case 1
                PS = 400e-9;
                primary_dead = 50e-9;
                secondary_dead = 13.33e-9;
                Power = (Ts/2)-PS-primary_dead-secondary_dead;
                R1 = 0.141;
                
            case 6
                PS = 250e-9;
                primary_dead = 110e-9;
                secondary_dead = 13.33e-9;
                Power = (Ts/2)-PS-primary_dead-secondary_dead;
                R1 = 0.212955;
                
            case 9
                PS = 410e-9;
                primary_dead = 86.66e-9;
                secondary_dead = 13.33e-9;
                Power = (Ts/2)-PS-primary_dead-secondary_dead;
                R1 = 0.16;
                
            case 12
                PS = 620e-9;
                primary_dead = 86.66e-9;
                secondary_dead = 13.33e-9;
                Power = (Ts/2)-PS-primary_dead-secondary_dead;
                R1 = 0.1277;
        end


%% There are several variable that must be filled out here:

% You must define u (Input variables in the order of [Independent
% Voltage Soruces MOSFET Forward Votlage (in order of netlist)
% Independent Current Sources]'
%%%% Example u  = [Vg Vfwd1 Vfwd2 Iout]';
 u = [Vg 1 1 1 1 1 1 1 1]';


% Define the inital guess of time intervals. This only defines the
% active switching time. Do not account for diode switching times.
%%%% For example a synchronous Buck converter would be: 
%%%%% ts = [Ts*(D-dead) Ts*(dead) Ts*(1-D-dead) Ts*(dead)];
%%%% But a non-synchronous buck covnerter would be
%%%%% ts = [Ts*(D-dead) Ts*(1-(D-dead))];

ts = [primary_dead PS secondary_dead Power primary_dead PS secondary_dead Power ];% The inital guess of time intervals



% List all of the numerical components in the netlist file for all
% FETs you must use the syntax used below:
          Numerical_Components = {'V1' Vg
        'C1' C1
        'C2' C2
        'L1' L1
        'L2' L2
        'L3' L3
        'M1_C' M1_C
        'M2_C' M2_C
        'M3_C' M3_C
        'M4_C' M4_C
        'M5_C' M5_C
        'M6_C' M6_C
        'M7_C' M7_C
        'M8_C' M8_C
        'R1' R1
        'R2' R2
        'R3' R3 
        'M1_R' 0.05
        'M2_R' 0.05
        'M3_R' 0.05
        'M4_R' 0.05
        'M5_R' 0.0015
        'M6_R' 0.0015
        'M7_R' 0.0015
        'M8_R' 0.0015
        
        };

% List out all char variables in the converter
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
Switch_Sequence = [
        OFF OFF OFF OFF ON OFF OFF ON %primary sw
        OFF ON ON OFF ON OFF OFF ON % phase shift
        OFF ON ON OFF OFF OFF OFF OFF % secondary sw
        OFF ON ON OFF OFF ON ON OFF % Reverse power
        OFF OFF OFF OFF OFF ON ON OFF % primary sw
        ON OFF OFF ON OFF ON ON OFF % phase shift
        ON OFF OFF ON OFF OFF OFF OFF % secondary sw
        ON OFF OFF ON ON OFF OFF ON]; % POWER


% List all the resistances of the diodes or FETS when they are on or
% off in order of Switch Variables as listed above
SW_OFF = ones(1,8).*10000000;

SW_ON = [0.05 0.05 0.05 0.05 0.0015 0.0015 0.0015 0.0015];

SW = [SW_OFF;SW_ON;SW_ON];

Diode_Forward_Voltage = [1 1 1 1 1 1 1 1]';



%% This creates and runs the parser to set up the converter

parse = NetListParse(); 
parse.initialize(filename,Voltage,Current);
parse.Component_Values = Numerical_Components; %%% Added to test numerical stuff in code
parse.read_file_num();
parse.
parse.ABCD_num();

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













