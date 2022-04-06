

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

X = [4 4 4 4 4 4 4 4 1.2315e6 1 1 1 1];


penalty  = 0;
%adjust = [ones(1,8), 1e-6, 100, 100, 100, 100] ;
%X = X./adjust;

FET_selection = X(1:8);
fs = X(end-4);
dt1 = X(end-3);
dt2 = X(end-2);
dt3 = X(end-1);
dt4 = X(end);
ron = [];
Coss = [];
tot_area = 0;

%% Determining Coss and ron based on w

for select_FET = 1:length(FET_selection)
    % Use typical value for Coss and max value for rds from data sheet
    % for inial look
    switch FET_selection(select_FET)
        
        %% EPC 2023
        case 1
            
            ron(select_FET) = 1.45e-3;
            Coss(select_FET) = 1530e-12;
            w = 6.05;
            l = 2.3;
            
            %% EPC 2014C
        case 2
            
            ron(select_FET) = 16e-3;
            Coss(select_FET) = 1530e-12;
            w = 1.7;
            l = 1.1;
            
            %% EPC 2015C
        case 3
            
            ron(select_FET) = 4e-3;
            Coss(select_FET) = 710e-12;
            w = 4.1;
            l = 1.6;
            
            %% EPC 2055
        case 4
            
            ron(select_FET) = 3.6e-3;
            Coss(select_FET) = 408e-12;
            w = 2.5;
            l = 1.5;
            
            %% EPC 2030
        case 5
            
            ron(select_FET) = 2.4e-3;
            Coss(select_FET) = 1120e-12;
            w = 4.6;
            l = 2.6;
            
            %% EPC 2024
        case 6
            
            ron(select_FET) = 1.5e-3;
            Coss(select_FET) = 1620e-12;
            w = 6.05;
            l = 2.3;
            
            %% EPC 2031
        case 7
            
            ron(select_FET) = 2.6e-3;
            Coss(select_FET) = 980e-12;
            w = 4.6;
            l = 2.6;
            
            
            %% EPC 2020
        case 8
            
            ron(select_FET) = 2.2e-3;
            Coss(select_FET) = 1020e-12;
            w = 6.05;
            l = 2.3;
            
            %% EPC 2065
        case 9
            
            ron(select_FET) = 3.6e-3;
            Coss(select_FET) = 534e-12;
            w = 3.5;
            l = 1.95;
            
            %% EPC 2029
        case 10
            
            ron(select_FET) = 3.2e-3;
            Coss(select_FET) = 820e-12;
            w = 4.6;
            l = 2.6;
            
            
            %% EPC 2021
        case 11
            
            ron(select_FET) = 2.2e-3;
            Coss(select_FET) = 1100e-12;
            w = 6.05;
            l = 2.3;
            
        otherwise
            
            stick = 100;
            return
            
    end
    tot_area = tot_area+w*l;
end

if 2*tot_area > (1.00001e-6)*100
    
    penalty = 2*tot_area*1e6;
    
else
    penalty = 0;
    
end


%% This is all caluclations to set up the variables need to find the SS
% Solution
Vg = 11.35;

Vin  = Vg;
Vb1 = 4.24;
Vb2 = 4.43;

VgPOS = 1;
Vb1POS = 2;
Vb2POS = 3;



Rb = 5e-3;
RL = 5e-3;


Co = 2e-6; ESRo = 2e-3;
Cfly1 = 5e-6; ESR1 = 3e-3;  % 5V Cap
Cfly2 = 2.5e-6; ESR2 = 3e-3; % 10V Cap
Lc = 100e-9;
Ts = 1/fs;

u = [Vg Vb1 Vb2 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5]';

modSchemes(:,:,1) = [
    0     1     1     1     0     1     0     1     1     0     0     1     1     1     0     1
    0     1     1     1     0     1     0     1     0     0     0     1     1     1     0     1
    0     1     1     1     0     1     0     1     0     1     0     1     1     1     0     1
    0     1     1     1     0     1     0     1     0     0     0     1     1     1     0     1
    0     1     1     1     0     1     0     1     1     0     0     1     1     1     0     1
    0     1     1     1     0     1     0     0     1     0     0     1     1     1     0     1
    0     1     1     1     0     1     1     0     1     0     0     1     1     1     0     1
    0     1     1     1     0     1     0     0     1     0     0     1     1     1     0     1
    ];

modSchemes(:,:,2) = [
    0     1     1     1     0     1     0     1     0     1     0     1     1     1     0     1
    0     1     1     1     0     1     0     0     0     1     0     1     1     1     0     1
    0     1     1     1     0     1     1     0     0     1     0     1     1     1     0     1
    0     1     1     1     0     1     1     0     0     0     0     1     1     1     0     1
    0     1     1     1     0     1     1     0     1     0     0     1     1     1     0     1
    0     1     1     1     0     1     1     0     0     0     0     1     1     1     0     1
    0     1     1     1     0     1     1     0     0     1     0     1     1     1     0     1
    0     1     1     1     0     1     0     0     0     1     0     1     1     1     0     1
    ];

modSchemes(:,:,3) = [
    0     1     1     1     0     1     1     0     0     1     0     0     1     0     0     1
    0     1     1     0     0     0     1     0     0     1     0     0     1     0     0     1
    0     1     1     0     1     0     1     0     0     1     0     0     1     0     0     1
    0     1     1     0     0     0     1     0     0     1     0     0     1     0     0     1
    0     1     1     1     0     1     1     0     0     1     0     0     1     0     0     1
    0     0     0     1     0     1     1     0     0     1     0     0     1     0     0     1
    1     0     0     1     0     1     1     0     0     1     0     0     1     0     0     1
    0     0     0     1     0     1     1     0     0     1     0     0     1     0     0     1
    ];

modSchemes(:,:,4) = [
    0     1     1     0     1     0     1     0     0     1     0     1     1     0     0     1
    0     1     1     0     0     0     1     0     0     1     0     1     1     0     0     0
    0     1     1     1     0     1     1     0     0     1     0     1     1     0     1     0
    0     0     0     1     0     1     1     0     0     1     0     0     1     0     0     0
    1     0     0     1     0     1     1     0     0     1     0     0     1     1     0     1
    0     0     0     1     0     1     1     0     0     1     0     0     0     1     0     1
    0     1     1     1     0     1     1     0     0     1     1     0     0     1     0     1
    0     1     1     0     0     0     1     0     0     1     0     0     0     0     0     1
    ];

modSchemes(:,:,5) = [
    0     1     1     1     0     1     1     0     0     1     0     1     1     0     1     0
    0     1     1     0     0     0     1     0     0     1     0     1     1     0     1     0
    0     1     1     0     1     0     1     0     0     1     0     1     1     0     1     0
    0     1     1     0     0     0     1     0     0     1     0     0     0     0     0     0
    0     1     1     1     0     1     1     0     0     1     1     0     0     1     0     1
    0     0     0     1     0     1     1     0     0     1     1     0     0     1     0     1
    1     0     0     1     0     1     1     0     0     1     1     0     0     1     0     1
    0     0     0     1     0     1     1     0     0     1     0     0     0     0     0     0
    ];

modSchemes(:,:,6) = [
    0     1     1     0     1     0     1     0     0     1     0     1     1     0     1     0
    0     1     1     0     0     0     1     0     0     1     0     0     0     0     1     0
    0     1     1     1     0     1     1     0     0     1     1     0     0     0     1     0
    0     0     0     1     0     1     1     0     0     1     1     0     0     0     0     0
    1     0     0     1     0     1     1     0     0     1     1     0     0     1     0     1
    0     0     0     1     0     1     1     0     0     1     1     0     0     0     0     0
    0     1     1     1     0     1     1     0     0     1     1     0     0     0     1     0
    0     1     1     0     0     0     1     0     0     1     0     0     0     0     1     0
    ];





%% FET Pairs
%{
M1 and M5    (1)
M2 and M4    (2)
M3 and M6    (3)
M9 and M12   (4)
M10 and M11  (5)
M13 and M17  (6)
M14 and M16  (7)
M15 and M18  (8)
%}


%% There are several variable that must be filled out here:

% You must define u (Input variables in the order of [Independent
% Voltage Soruces MOSFET Forward Votlage (in order of netlist)
% Independent Current Sources]'
%%%% Example u  = [Vg Vfwd1 Vfwd2 Iout]';
% Assigned later


% Define the inital guess of time intervals. This only defines the
% active switching time. Do not account for diode switching times.
%%%% For example a synchronous Buck converter would be:
%%%%% ts = [Ts*(D-dead) Ts*(dead) Ts*(1-D-dead) Ts*(dead)];
%%%% But a non-synchronous buck covnerter would be
%%%%% ts = [Ts*(D-dead) Ts*(1-(D-dead))];
d = 0.5;
dt1 = 0.00370/2;
dt2 = 0.00370/2;
dt3 = 0.00370/2;
dt4 = 0.00370/2;
ds = [d-dt1, dt1, 1-d-dt2, dt2, d-dt3, dt3, 1-d-dt4, dt4];
ts = ds/sum(ds)*Ts;


% The inital guess of time intervals % The inital guess of time intervals
% Assigned later dynamically


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
    'M1_C' Coss(1)
    'M2_C' Coss(2)
    'M3_C' Coss(3)
    'M4_C' Coss(2)
    'M5_C' Coss(1)
    'M6_C' Coss(3)
    'M7_C' Coss(4)
    'M8_C' Coss(5)
    'M9_C' Coss(5)
    'M10_C' Coss(4)
    'M11_C' Coss(6)
    'M12_C' Coss(7)
    'M13_C' Coss(8)
    'M14_C' Coss(7)
    'M15_C' Coss(6)
    'M16_C' Coss(8)
    'M1_R' ron(1)
    'M2_R' ron(2)
    'M3_R' ron(3)
    'M4_R' ron(2)
    'M5_R' ron(1)
    'M6_R' ron(3)
    'M7_R' ron(4)
    'M8_R' ron(5)
    'M9_R' ron(5)
    'M10_R' ron(4)
    'M11_R' ron(6)
    'M12_R' ron(7)
    'M13_R' ron(8)
    'M14_R' ron(7)
    'M15_R' ron(6)
    'M16_R' ron(8)
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
SW_ON = [ron(1) ron(2) ron(3) ron(2) ron(1) ron(3) ron(4) ron(5) ron(5) ron(4) ron(6) ron(7) ron(8) ron(7) ron(6) ron(8)];
SW = [SW_OFF;SW_ON;SW_ON];

Diode_Forward_Voltage = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]'.*1.5;




%% This creates and runs the parser to set up the converter

parse = NetListParse();
parse.initialize(filename,Voltage,Current);
parse.Component_Values = Numerical_Components; %%% Added to test numerical stuff in code
parse.cutset_loop_num();

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
sim = SMPSim;



etaSim = [];
PlossSim = [];
PoutSim = [];
conditions = [];
PossSim = [];
drange = linspace(0.05, .95, 8);
Vscloc1 = 38;
Vscloc2 = 39;
Illoc = 1;
Ib1loc = 2;
Ib2loc = 3;

sim=Run_SS_Converter_num_no_diode(sim,conv);

for k = 1:1:length(parse.StateNumbers)
    if strcmp(parse.OutputNamesCD{parse.StateNumbers(k),1}(1,end),'A')
        top.stateLabels(end+1,1) = strcat('I_{', parse.StateNames(k,1),'} (A)');
        top.stateLabels_Opp(end+1,1) = strcat('V_{', parse.StateNames(k,1),'} (V)');
    else
        top.stateLabels(end+1,1) = strcat(parse.OutputNamesCD{parse.StateNumbers(k),1}(1,end),'_{', parse.StateNames(k,1),'} (V)');
        top.stateLabels_Opp(end+1,1) = strcat('I_{', parse.StateNames(k,1),'} (A)');
    end
end

conv.As = sim.As;
conv.Bs = sim.Bs;
conv.Cs = sim.Cs;
conv.Ds = sim.Ds;
conv.eigA = sim.eigA;


sim.As = conv.As;
sim.Bs = conv.Bs;
sim.Cs = conv.Cs;
sim.Ds = conv.Ds;
sim.eigA = conv.eigA;
sim.ts = conv.ts;

Xss = sim.SS_Soln();


parse.find_diode_new(conv.order,conv.Switch_Sequence,conv.Fwd_Voltage)
iterations = 50;
cycle = 0;
fail = 1;
%fprintf('--------------\n')
fail = sim.Three_tier_diode_correct_num(iterations,0,0);

while fail && cycle < 2
    fail = sim.Three_tier_diode_correct_num(iterations,0,1);
    cycle = cycle+1;
end


Xss = sim.SS_Soln();
[ avgXs, avgYs ] = sim.ssAvgs(Xss);

Poss = 0;

% for when there is no rms calculation
%%{
eta = (avgYs(Ib1loc)*Vb1 + avgYs(Ib2loc)*Vb2 - Poss)/(Vin*-avgYs(Illoc) - avgYs(Ib1loc)^2*Rb - avgYs(Ib2loc)^2*Rb);
Ploss = -(avgYs(Ib1loc)*Vb1 + avgYs(Ib2loc)*Vb2) + -(Vin*avgYs(Illoc)) + Poss;
Pout = (avgYs(Ib1loc)*Vb1 + avgYs(Ib2loc)*Vb2) - Poss;

toc

return



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













