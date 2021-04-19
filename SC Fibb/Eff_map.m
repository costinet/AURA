function [stick] = Eff_map(X)

% Inside GA Loop
% This is going to make a eff plot for the converter Init


%     _   _   _  ____    _
%    / \ | | | |/ _  |  / \
%   / _ \| | | | (_| | / _ \
%  / ___ | |_| |> _  |/ ___ \
% /_/   \_\___//_/ |_/_/   \_\

X = X./[1 1 1 1 1 1 1 1 20 20 20 20 20 20 20 20 1e-6];
FET_selection = X(1:8);
W_selection = X(9:16);
fs = X(17);

%% Scaling
% Frequency / 1e6
% Width *20
% FET Selection same


% Current konwn issues use FETs for diodes. (Just have them never turn
% on)


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

ron = [];
Coss = [];
tot_area = 0;
%% Determining Coss and ron based on w

for select_FET = 1:8
    
    switch FET_selection(select_FET)
        
        %% 8HVnLDMOS nbl
        case 1
            a =  0.001378;
            b = -1;
            ron(select_FET) = a*W_selection(select_FET)^b;
            
            p1 = 2.925e-9;
            p2 = -2.787e-12;
            Coss(select_FET) = p1*W_selection(select_FET)+p2;
            
            L1=800*10^-9;
            
            tot_area = tot_area+W_selection(select_FET)*L1;
            
            
            %% 8HVnLDMOS iso
        case 2
            a =  0.001453;
            b = -1;
            ron(select_FET) = a*W_selection(select_FET)^b;
            
            p1 = 2.093e-9;
            p2 = -1.058e-12;
            Coss(select_FET) = p1*W_selection(select_FET)+p2;
            
            L2=900*10^-9;
            
            tot_area = tot_area+W_selection(select_FET)*L2;
            
            
            %% 12HVnLDMOS iso hp mac
        case 3
            a =  0.001447;
            b = -1;
            ron(select_FET) = a*W_selection(select_FET)^b;
            
            p1 = 2.487e-9;
            p2 = -6.393e-13;
            Coss(select_FET) = p1*W_selection(select_FET)+p2;
            
            L3=900*10^-9;
            
            tot_area = tot_area+W_selection(select_FET)*L3;
            
            
            %% 12HVnLDMOS iso mac
        case 4
            a =  0.001656;
            b = -1;
            ron(select_FET) = a*W_selection(select_FET)^b;
            
            p1 = 2.281e-9;
            p2 = -3.1e-12;
            Coss(select_FET) = p1*W_selection(select_FET)+p2;
            
            L4=1*10^-6;
            
            
            tot_area = tot_area+W_selection(select_FET)*L4;
            
            %% 12HVnLDMOS nbl hp mac
        case 5
            a =  0.001447;
            b = -1;
            ron(select_FET) = a*W_selection(select_FET)^b;
            
            p1 = 2.487e-9;
            p2 = -6.393e-13;
            Coss(select_FET) = p1*W_selection(select_FET)+p2;
            
            L5=0.9*10^-6;
            
            tot_area = tot_area+W_selection(select_FET)*L5;
            
            
            %% 12HVnLDMOS nbl mac
        case 6
            a =  0.001585;
            b = -1;
            ron(select_FET) = a*W_selection(select_FET)^b;
            
            p1 = 3.176e-9;
            p2 = -7.226e-12;
            Coss(select_FET) = p1*W_selection(select_FET)+p2;
            
            L6=1*10^-6;
            
            tot_area = tot_area+W_selection(select_FET)*L6;
            
            %% 12HVnLDMOS nbl mr mac
        case 7
            a =  0.002556;
            b = -1;
            ron(select_FET) = a*W_selection(select_FET)^b;
            
            p1 = 2.828e-9;
            p2 = -1.703e-12;
            Coss(select_FET) = p1*W_selection(select_FET)+p2;
            
            L7=.9*10^-6;
            
            tot_area = tot_area+W_selection(select_FET)*L7;
            
    end
    
end


%% FET Pairs
%{
M1 and M5
M2 and M4
M3 and M6
M7 and M10
M8 and M9
M11 and M15
M12 and M14
M13 and M16
%}

[Coss1, Coss5]=deal(Coss(1));
[Coss2, Coss4]=deal(Coss(2));
[Coss3, Coss6]=deal(Coss(3));
[Coss7, Coss10]=deal(Coss(4));
[Coss8, Coss9]=deal(Coss(5));
[Coss11, Coss15]=deal(Coss(6));
[Coss12, Coss14]=deal(Coss(7));
[Coss13, Coss16]=deal(Coss(8));

[ron1, ron5]=deal(ron(1));
[ron2, ron4]=deal(ron(2));
[ron3, ron6]=deal(ron(3));
[ron7, ron10]=deal(ron(4));
[ron8, ron9]=deal(ron(5));
[ron11, ron15]=deal(ron(6));
[ron12, ron14]=deal(ron(7));
[ron13, ron16]=deal(ron(8));



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

%fs = 2e6;
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
    'M1_C' Coss1
    'M2_C' Coss2
    'M3_C' Coss3
    'M4_C' Coss4
    'M5_C' Coss5
    'M6_C' Coss6
    'M7_C' Coss7
    'M8_C' Coss8
    'M9_C' Coss9
    'M10_C' Coss10
    'M11_C' Coss11
    'M12_C' Coss12
    'M13_C' Coss13
    'M14_C' Coss14
    'M15_C' Coss15
    'M16_C' Coss16
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
SW_ON = [ron1 ron2 ron3 ron4 ron5 ron6 ron7 ron8 ron9 ron10 ron11 ron12 ron13 ron14 ron15 ron16];
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

for i = 1:1:length(parse.StateNumbers)
    if strcmp(parse.OutputNamesCD{parse.StateNumbers(i),1}(1,end),'A')
        top.stateLabels(end+1,1) = strcat('I_{', parse.StateNames(i,1),'} (A)');
        top.stateLabels_Opp(end+1,1) = strcat('V_{', parse.StateNames(i,1),'} (V)');
    else
        top.stateLabels(end+1,1) = strcat(parse.OutputNamesCD{parse.StateNumbers(i),1}(1,end),'_{', parse.StateNames(i,1),'} (V)');
        top.stateLabels_Opp(end+1,1) = strcat('I_{', parse.StateNames(i,1),'} (A)');
    end
end

parse.find_diode_new(Order,conv.Switch_Sequence,Diode_Forward_Voltage);


Vgloc = 1;

simulator = SMPSim;
etaSim = [];
PlossSim = [];
PoutSim = [];
conditions = [];
PossSim = [];

drange = linspace(0.01, 0.99, 10);

h = waitbar(0, 'Simulating');

for i = 1:size(modSchemes,3)
    conv.Switch_Sequence = modSchemes(:,:,i);
    parse.find_diode_new(Order,conv.Switch_Sequence,Diode_Forward_Voltage);
    waitbar((i-1)/size(modSchemes,3), h);
    for d = drange
        ds = [d, 1-d, d, 1-d];
        ts = ds/sum(ds)*Ts;
        conv.ts = ts;
        
        %% Find zero-power Vin
        conv.Element_Properties(find(strcmp(conv.Element_Properties, 'R3')),2) = {1e6};
        EVAL_SS_Converter_Fib(conv)
        simulator.loadTestConverter2(conv);
        simulator.eigA = parse.eigA;
        simulator.binary = conv.Switch_Sequence;
        Xss = simulator.SS_Soln();

        [ avgX, avgY ] = simulator.ssAvgs(Xss);
        
        OLVin = avgY(38)+avgY(39);
        MaxVin = OLVin+1;
        
        if OLVin<5
            continue
        end
        
        
        
        conv.Element_Properties(find(strcmp(conv.Element_Properties, 'R3')),2) = {Rb};
        EVAL_SS_Converter_Fib(conv)
        simulator.loadTestConverter2(conv);
        simulator.eigA = parse.eigA;
        simulator.binary = conv.Switch_Sequence;
        Xss = simulator.SS_Soln();
        
         waitbar((i-1)/size(modSchemes,3) + 1/size(modSchemes,3)*(find(drange==d)-1)/length(drange) , h);
        
        
        
        Vinrange = OLVin:.02:MaxVin;
        for Vin = Vinrange
            simulator.u(Vgloc) = Vin;
            
            Xss = simulator.SS_Soln();
            [ avgX, avgY ] = simulator.ssAvgs(simulator.Xs);
            eff=-((avgY(3)*avgY(30))+(avgY(2)*avgY(29)))/(avgY(1)*avgY(28));
            Pout = ((avgY(3)*avgY(30))+(avgY(2)*avgY(29)));
            Ploss = (Pout/eff)-Pout;
            
            if eff>1||eff<.8
                continue
            end
            
            if Pout>100|| Pout<0
                continue
            end
            
            PlossSim = [PlossSim Ploss];
            etaSim = [etaSim; eff];
            PoutSim = [PoutSim; Pout];
            conditions = [conditions; d, i, Vin];
            
            
            
        end
    end
end

close(h)


J = 45645314565;


mineff = .5;
etaSim(etaSim<mineff) = mineff;
etaSim(etaSim>1) = mineff;

locs = etaSim > mineff;

VgSim = conditions(locs,end);
F = scatteredInterpolant(VgSim, PoutSim(locs), etaSim(locs)*100, 'linear','none');

figure(1)

Vgrange = 5:.25:22;
PoutRange = 0:1:80;
[VgMesh, PoutMesh] = meshgrid(Vgrange, PoutRange);
[f,c] = contourf(Vgrange,PoutRange,F(VgMesh,PoutMesh),'ShowText','on');
xlabel('Vg');
ylabel('P_{out}');
% ylim([0 80])
colorbar
% c.LevelList = [.8 .9 .95 .96 .97 .98 .99]*100;

hold on;

figure(2)


% 
% figure(2)
% subplot(2,1,1)
% F = scatteredInterpolant(VgSim, PoutSim(locs), PlossSim(locs), 'linear','none');
% [f,c] = contourf(Vgrange,PoutRange,F(VgMesh,PoutMesh),'ShowText','on');
% levels = c.LevelList;
% colorbar
% 
% subplot(2,1,2)
% F = scatteredInterpolant(VgSim, PoutSim(locs), PossSim(locs), 'linear','none');
% [f,c] = contourf(Vgrange,PoutRange,F(VgMesh,PoutMesh),'ShowText','on');
% colorbar
% 
F2 = scatteredInterpolant(VgSim, PoutSim(locs), conditions(locs,2), 'linear','none');
contourf(Vgrange,PoutRange,F2(VgMesh,PoutMesh),'ShowText','on');
scatter(VgSim, PoutSim(locs),[],conditions(locs,2));

[c] = scatter(conditions(:,end),PoutSim,[],etaSim*100);
xlabel('Vg');
ylabel('P_{out}');
ylim([0 100])
colorbar
%}
Added_eff = max(etaSim(conditions(:,3)>18 & PoutSim<20)) + ... 
   max(etaSim(conditions(:,3)>18 & PoutSim>60)) + ...
   max(etaSim(conditions(:,3)>18 & PoutSim<60 & PoutSim>40)) + ...
   max(etaSim(conditions(:,3)>18 & PoutSim<40 & PoutSim>20)) + ...
   max(etaSim(conditions(:,3)<10 & PoutSim<20)) + ... 
   max(etaSim(conditions(:,3)<10 & PoutSim>60)) + ...
   max(etaSim(conditions(:,3)<10 & PoutSim<60 & PoutSim>40)) + ...
   max(etaSim(conditions(:,3)<10 & PoutSim<40 & PoutSim>20)) + ...
   max(etaSim(conditions(:,3)>10 & conditions(:,3)<14 & PoutSim<20)) + ... 
   max(etaSim(conditions(:,3)>10 & conditions(:,3)<14 & PoutSim>60)) + ...
   max(etaSim(conditions(:,3)>10 & conditions(:,3)<14 & PoutSim<60 & PoutSim>40)) + ...
   max(etaSim(conditions(:,3)>10 & conditions(:,3)<14 & PoutSim<40 & PoutSim>20)) + ...
   max(etaSim(conditions(:,3)>14 & conditions(:,3)<18 & PoutSim<20)) + ... 
   max(etaSim(conditions(:,3)>14 & conditions(:,3)<18 & PoutSim>60)) + ...
   max(etaSim(conditions(:,3)>14 & conditions(:,3)<18 & PoutSim<60 & PoutSim>40)) + ...
   max(etaSim(conditions(:,3)>14 & conditions(:,3)<18 & PoutSim<40 & PoutSim>20));

stick = Added_eff/16;
stick = (1-stick)*100;
% That's all Folk's
end            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            