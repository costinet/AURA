
%     _   _   _  ____    _    
%    / \ | | | |/ _  |  / \   
%   / _ \| | | | (_| | / _ \  
%  / ___ | |_| |> _  |/ ___ \ 
% /_/   \_\___//_/ |_/_/   \_\



% Current konwn issues use FETs for diodes. (Just have them never turn
% on)

tic

% Place in the filename
filename = 'Buck_Boost.net'; % Place netlist filename here that you want to run

% Can also place component values here to get their voltage and
% current output waveforms such as: 

Voltage = {'V1'
    };
Current = {'V1'
    };


%% This is all caluclations to set up the variables need to find the SS
% Solution




        Vg = 5;
        Pout = 25;
        Dbuck = 0.92;
        Dboost = 0.3118;
        M = Dbuck/(1-Dboost);
        Vout = M*Vg;
        Iout = Pout/Vout;
        R1 = .73+.73;
        L1 = 200e-9;
        C1 = 23.5e-6;
        C2 = 4.4e-6;
        
        % EPC 2024
        r_on = 1.5e-3;
        Cds = 1620e-12;
        [M1_R,M2_R,M3_R,M4_R,M5_R,M6_R] = deal(r_on);
        [M1_C,M2_C,M3_C,M4_C,M5_C,M6_C] = deal(Cds);
        
        
        Numerical_Components = {'V1' Vg
            'C1' C1
            'C2' C2
            'L1' L1
            'M1_C' M1_C
            'M2_C' M2_C
            'M3_C' M3_C
            'M4_C' M4_C
            'M5_C' M5_C
            'M6_C' M6_C
            'R1' R1 };
        
        
        Switch_Resistors = {'M1_R'
            'M2_R'
            'M3_R'
            'M4_R'
            'M5_R'
            'M6_R'};
        
        fs = 1000e3;
        Ts = 1/fs;
        dt = 5e-9;
        
        ON = 1;
        OFF = 0;
        
        
        SW_OFF = ones(1,6).*10000000;
        
        SW_ON = [r_on r_on r_on r_on r_on r_on];
        
        SW = [SW_OFF;SW_ON;SW_ON];
        
        Order = [1 2 3 4 5 6 7 8 9 10 11 12];
        
        
        u = [Vg 1 1 1 1 1 1]';
        
        
        % Duty cycle is in pu

        
        
        ts = [Dboost*Ts-dt  dt  (Dbuck/2-Dboost)*Ts-dt  dt  (0.5-Dbuck/2)*Ts-dt  dt  Dboost*Ts-dt  dt  (Dbuck/2-Dboost)*Ts-dt  dt  (0.5-Dbuck/2)*Ts-dt  dt];
        
        Binary_for_DAB = [
         ON  OFF OFF ON  OFF ON   % Q6 ON
         ON  OFF OFF ON  OFF OFF  % Dead Q6 Q3
         ON  OFF ON  ON  OFF OFF  % Q3 ON
         OFF OFF ON  ON  OFF OFF  % Dead Q1 Q2
         OFF ON  ON  ON  OFF OFF  % Q2 ON
         OFF OFF ON  OFF OFF OFF  % Dead
         
         ON  OFF ON  OFF ON  OFF  % Q5 Q1 ON
         ON  OFF ON  OFF OFF OFF  % Dead Q5 Q4
         ON  OFF ON  ON  OFF OFF  % Q4 ON
         OFF OFF ON  ON  OFF OFF  % Dead Q1 Q2         
         OFF ON  ON  ON  OFF OFF  % Q2 ON
         OFF OFF OFF ON  OFF OFF  % Dead
        ];
         
    
M1_C_Resist = 10000;


%% There are several variable that must be filled out here:

% You must define u (Input variables in the order of [Independent
% Voltage Soruces MOSFET Forward Votlage (in order of netlist)
% Independent Current Sources]'
%%%% Example u  = [Vg Vfwd1 Vfwd2 Iout]';
u = [Vg 1 1 1 1 1 1]';


% Define the inital guess of time intervals. This only defines the
% active switching time. Do not account for diode switching times.
%%%% For example a synchronous Buck converter would be: 
%%%%% ts = [Ts*(D-dead) Ts*(dead) Ts*(1-D-dead) Ts*(dead)];
%%%% But a non-synchronous buck covnerter would be
%%%%% ts = [Ts*(D-dead) Ts*(1-(D-dead))];

 ts = ts; % The inital guess of time intervals % The inital guess of time intervals



% List all of the numerical components in the netlist file for all
% FETs you must use the syntax used below:
        Numerical_Components = {'C1' C1
            'C2' C2
            'L1' L1
            'M1_C' M1_C
            'M2_C' M2_C
            'M3_C' M3_C
            'M4_C' M4_C
            'M5_C' M5_C
            'M6_C' M6_C
            'R1' R1

            };

% List out all char variables in the 
       Switch_Resistors = {'M1_R'
            'M2_R'
            'M3_R'
            'M4_R'
            'M5_R'
            'M6_R'};
% List of the switch sequency. Organized by: the FETs (column) vs time
% interval (rows) matching Switch_Resistors and ts respectivly

ON = 1;
OFF = 0;
Switch_Sequence = [
    ON  OFF OFF ON  OFF ON   % Q6 ON
    ON  OFF OFF ON  OFF OFF  % Dead Q6 Q3
    ON  OFF ON  ON  OFF OFF  % Q3 ON
    OFF OFF ON  ON  OFF OFF  % Dead Q1 Q2
    OFF ON  ON  ON  OFF OFF  % Q2 ON
    OFF OFF ON  OFF OFF OFF  % Dead
    
    ON  OFF ON  OFF ON  OFF  % Q5 Q1 ON
    ON  OFF ON  OFF OFF OFF  % Dead Q5 Q4
    ON  OFF ON  ON  OFF OFF  % Q4 ON
    OFF OFF ON  ON  OFF OFF  % Dead Q1 Q2
    OFF ON  ON  ON  OFF OFF  % Q2 ON
    OFF OFF OFF ON  OFF OFF  % Dead
    ];


% List all the resistances of the diodes or FETS when they are on or
% off
SW_OFF = ones(1,6).*10000000;

SW_ON = [M1_R,M2_R,M3_R,M4_R,M5_R,M6_R];

SW = [SW_OFF;SW_ON;SW_ON];


Diode_Forward_Voltage = [1 1 1 1 1 1]';



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





         