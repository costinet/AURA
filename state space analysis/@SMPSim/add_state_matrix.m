function [A,B,C,D,eigA] = add_state_matrix(obj,new_state)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


%% DAB Converter
%{
Vg = 48;
L3 = 1e-6;
L2 = 1296e-6;
L1 = 76e-6;
C1 = (4*100+680)*10^-6; %Cout
% C1 = (4*100)*10^-6; %Cout
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




M1_R_D = 0.05; % ronHS
M2_R_D = 0.05; % ronLS
M3_R_D = 0.05; % ronHS
M4_R_D = 0.05; % ronLS
M5_R_D = .0015; % ronHS
M6_R_D = .0015; % ronLSR2 = 0.16;
R3 = 0.001;

%R1 =  0.16; % the output resistor

%R1 = 0.05; % Output resistor of test
 % Output resistor of test

R2 = 4.57535; % the inductor resistance % 1 ohm on the primary 1Ohm on the secondary
dt = Ts/1000;%5e-10;
Vdr = 5;
M1_R_ON = 0.05; % ronHS
M2_R_ON = 0.05; % ronLS
M3_R_ON = 0.05; % ronHS
M4_R_ON = 0.05; % ronLS
M5_R_ON = 0.0015; % ronHS
M6_R_ON = .0015; % ronLS
M7_R_ON = .0015; % ronHS
M8_R_ON = .0015; % ronLS
M7_R_D = .0015; % ronHS
M8_R_D = .0015; % ronLS



[D1_R_OFF,M1_R_OFF,M2_R_OFF,M3_R_OFF,M4_R_OFF,M5_R_OFF,M6_R_OFF,M7_R_OFF,M8_R_OFF] = deal(100000000); % ronLS

Order = [1 2 3 4 5 6 7 8]; % The order that the states must go in after being parsed
% ts = [Ts*.2376 Ts*.0277 Ts*.23 Ts*.00023 Ts*.2376 Ts*.02775 Ts*.23 Ts*.00023];



%[D1_R,M1_R,M2_R,M3_R,M4_R,M5_R,M6_R,M7_R,M8_R] = deal(100000000); % ronLS



% ts = [Ts*.33 Ts*.01 Ts*.16 Ts*.001 Ts*.33 Ts*.01 Ts*.16 Ts*.001];
chooses = 9;
switch chooses
    
    case 1
        PS = 400e-9;
        primary_dead = 50e-9;
        secondary_dead = 13.33e-9;
        Power = (Ts/2)-PS-primary_dead-secondary_dead;
        R1 = 0.141;
        
        case 6
        PS = 250-9;
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

ts = [primary_dead PS secondary_dead Power primary_dead PS secondary_dead Power ];

%%%%%%ts = [(Ts/2)-50e-9-400e-9-13.33e-9 50e-9 400e-9 13.33e-9 (Ts/2)-50e-9-400e-9-13.33e-9 50e-9 400e-9 13.33e-9];

%%%%%%ts = [50e-9 400e-9 13.33e-9 (Ts/2)-50e-9-400e-9-13.33e-9 50e-9 400e-9 13.33e-9 (Ts/2)-50e-9-400e-9-13.33e-9];



%ts = [2.4960e-06 100e-9 200e-9 100e-9 2.4960e-06 100e-9 200e-9 100e-9];


%ts = [2.1e-06 100e-9 200e-9 10e-11 2.1e-06 100e-9 200e-9 10e-11];

% ts = [2.1e-06 50-9 400e-9 10e-9 2.1e-06 50e-9 400e-9 10e-9];


u = [Vg 1 1 1 1 1 1 1 1]';


TestparseWaveform = false;

ON = 1;
OFF = 0;


SW_OFF = ones(1,8).*10000000;

SW_ON = [0.05 0.05 0.05 0.05 0.0015 0.0015 0.0015 0.0015];

SW = [SW_OFF;SW_ON;SW_ON];


%Out = [M1_R,M2_R,M3_R,M4_R,M5_R,M6_R,M7_R,M8_R];



% Binary_for_DAB = [
%     OFF ON ON OFF OFF ON ON OFF % Reverse power
%     OFF OFF OFF OFF OFF ON ON OFF % primary sw
%     ON OFF OFF ON OFF ON ON OFF % phase shift
%     ON OFF OFF ON OFF OFF OFF OFF % secondary sw
%     ON OFF OFF ON ON OFF OFF ON % POWER
%     OFF OFF OFF OFF ON OFF OFF ON %primary sw
%     OFF ON ON OFF ON OFF OFF ON % phase shift
%     OFF ON ON OFF OFF OFF OFF OFF]; % secondary sw

Binary_for_DAB = [

    OFF OFF OFF OFF OFF ON ON OFF % primary sw
    ON OFF OFF ON OFF ON ON OFF % phase shift
    ON OFF OFF ON OFF OFF OFF OFF % secondary sw
    ON OFF OFF ON ON OFF OFF ON % POWER
    OFF OFF OFF OFF ON OFF OFF ON %primary sw
    OFF ON ON OFF ON OFF OFF ON % phase shift
    OFF ON ON OFF OFF OFF OFF OFF % secondary sw
    OFF ON ON OFF OFF ON ON OFF]; % Reverse power

%}


%% New Buck Converter
%{
Vg = 12;
L1 = 1e-6; %L
C1 = (4*100+680)*10^-6; %Cout
fs = 0.2e6;
Ts = 1/fs;
V = 5;
%Io = 1; % was 1
M1_C = 1620e-12; % CHS
M2_C = 1620e-12; % LHS
R3 = 0.001;
R2 = 0.16;
R1 =  1; % Rl also known as R1 (one or L (lowercase))
dt = Ts/100;%5e-10;
Vdr = 5;
M1_R_ON = 0.0015; % ronHS
M2_R_ON = 0.0015; % ronLS
D1_R_D = 0.0015; % ronLS
M1_R_D = 0.0015; % ronHS
M2_R_D = 0.0015; % ronLS
[D1_R_OFF,M1_R_OFF,M2_R_OFF] = deal(100000000); % ronLS

%Order = [2 4]; % The order that the states must go in after being parsed
%ts = [Ts*.5 Ts*.5]; % The inital guess of time intervals


Order = [1 2 3 4]; % The order that the states must go in after being parsed
% ts = [Ts*.5-dt dt Ts*.5-dt dt]; % The inital guess of time intervals
% u = [Vg  1 1 Io]';


% Order = [1 2 3 4 5 6 7 8]; % The order that the states must go in after being parsed
% ts = [Ts*.2376 Ts*.0277 Ts*.23 Ts*.00023 Ts*.2376 Ts*.02775 Ts*.23 Ts*.00023];



%[D1_R,M1_R,M2_R,M3_R,M4_R,M5_R,M6_R,M7_R,M8_R] = deal(100000000); % ronLS

D = 5/12;
dead = 0.001/fs;

ts = [Ts*D-dead dead Ts*(1-D)-dead dead];

u = [Vg 1 1]';


ON = 1;
OFF = 0;


SW_OFF = ones(1,2).*10000000;

SW_ON = [M1_R_ON M2_R_ON];

SW = [SW_OFF;SW_ON;SW_ON];


%Out = [M1_R,M2_R,M3_R,M4_R,M5_R,M6_R,M7_R,M8_R];



Binary_for_DAB = [
    ON OFF
    OFF OFF
    OFF ON
    OFF OFF];

%}

%%  Hybrid Dixson Converter

%{
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

SW = [SW_OFF;SW_ON;SW_ON];

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
%}

%% MR-Buck Converter
%{
Vg = 12;
L1 = 1e-6; %L
C1 = 1*10^-6; %Cout
L2 = 20e-9; % Resonate inductor
fs = 10e6;
Ts = 1/fs;
V = 3;
Io = 1; % was 1
M1_C = 0.7e-9; % CHS
M2_C = 0.7e-9; % LHS
dt = Ts/100;%5e-10;
Vdr = 5;
M1_R_ON = 0.002; % ronHS
M2_R_ON = 0.002; % ronLS
D1_R_D = 0.002; % ronLS
M1_R_D = 0.002; % ronHS
M2_R_D = 0.002; % ronLS
[D1_R_OFF,M1_R_OFF,M2_R_OFF] = deal(100000000); % ronLS

%Order = [2 4]; % The order that the states must go in after being parsed
%ts = [Ts*.5 Ts*.5]; % The inital guess of time intervals


Order = [1 2 3 4]; % The order that the states must go in after being parsed
% ts = [Ts*.5-dt dt Ts*.5-dt dt]; % The inital guess of time intervals
% u = [Vg  1 1 Io]';


% Order = [1 2 3 4 5 6 7 8]; % The order that the states must go in after being parsed
% ts = [Ts*.2376 Ts*.0277 Ts*.23 Ts*.00023 Ts*.2376 Ts*.02775 Ts*.23 Ts*.00023];



%[D1_R,M1_R,M2_R,M3_R,M4_R,M5_R,M6_R,M7_R,M8_R] = deal(100000000); % ronLS

D = 5/12;
%dead = 0.001/fs;
dead = 0.1/fs;
ts = [25e-9 25e-9 25e-9 25e-9];

u = [Vg Io 1 1]';


TestparseWaveform = false;

ON = 1;
OFF = 0;


SW_OFF = ones(1,2).*10000000;

SW_ON = [M1_R_ON M2_R_ON];

SW = [SW_OFF;SW_ON;SW_ON];




%}



%% Dickson
%{
Vg = 48;
Iload = 5;
fs = 1e6;
Ts = 1/fs;
dt = Ts/50;
L1 = 500e-9;
C8 = 10e-6;
r_on = 4e-3;
Cds = 710e-12;
Cx = 2e-6;
Cx_ESR = 4e-3;
[M1_R,M2_R,M3_R,M4_R,M5_R,M6_R,M7_R,M8_R,M9_R,M10_R,M11_R,M12_R] = deal(r_on);
[M1_C,M2_C,M3_C,M4_C,M5_C,M6_C,M7_C,M8_C,M9_C,M10_C,M11_C,M12_C] = deal(Cds);
[C1,C2,C3,C4,C5,C6,C7] = deal(Cx);
[R1,R2,R3,R4,R5,R6,R7,R8] = deal(Cx_ESR);


% Order = [1 2]; % The order that the states must go in after being parsed
%  ts = [Ts*.5 Ts*.5]; % The inital guess of time intervals
u = [Vg  1 1 1 1 1 1 1 1 1 1 1 1 Iload]';
Order  = [1 2 3 4];
ts = [Ts*.5-dt dt Ts*.5-dt dt]; % The inital guess of time intervals
% u = [Vg  1 1 Io]';


ON = 1;
OFF = 0;


SW_OFF = ones(1,12).*10000000;

SW_ON = [M1_R,M2_R,M3_R,M4_R,M5_R,M6_R,M7_R,M8_R,M9_R,M10_R,M11_R,M12_R];

SW = [SW_OFF;SW_ON;SW_ON];

%}

%% Active Clamp Flyback
%{
Vg = 5;
Iload = 5;
fs = 0.5e6;
Ts = 1/fs;
dt = Ts/50;
L1 = 1e-3; % Primary turns
L2 = 16e-3; % Secondary turns
L3 = 5e-3;  % Magnetizing inductance
r_on = 16e-3;
Cds = 150e-12;
L4 = 5e-6; % Resonate indcutor
C2 = 5e-9; % Resonate Cap
C1 = 50e-6; % Output Cap
R1 = 133.3333;
[M1_R,M2_R,M3_R] = deal(r_on);
[M1_C,M2_C,M3_C] = deal(Cds);


% Order = [1 2]; % The order that the states must go in after being parsed
%  ts = [Ts*.5 Ts*.5]; % The inital guess of time intervals
 u = [Vg 1 1 1]';
Order  = [1 2 3 4];
 ts = [Ts*.5-dt dt Ts*.5-dt dt]; % The inital guess of time intervals
% u = [Vg  1 1 Io]';


ON = 1;
OFF = 0;

 Binary_for_DAB = [
ON OFF OFF
OFF OFF OFF
OFF ON OFF
OFF OFF OFF
     ];

SW_OFF = ones(1,3).*10000000;

SW_ON = [M1_R,M2_R,M3_R];

SW = [SW_OFF;SW_ON;SW_ON];

%}

%%  Flyback Converters
%{


Vg = 5;
Iload = 5;
fs = 0.5e6;
Ts = 1/fs;
dt = Ts/50;
L1 = 1e-3; % Primary turns
L2 = 16e-3; % Secondary turns
L3 = 5e-3;  % Magnetizing inductance
r_on = 16e-3;
Cds = 150e-12;
L4 = 5e-6; % Resonate indcutor
C2 = 5e-9; % Resonate Cap
C1 = 50e-6; % Output Cap
R1 = 133.3333;
[M1_R,M2_R] = deal(r_on);
[M1_C,M2_C] = deal(Cds);

% Order = [1 2]; % The order that the states must go in after being parsed
%  ts = [Ts*.5 Ts*.5]; % The inital guess of time intervals
 u = [Vg 1 1 ]';
Order  = [1 2];
 ts = [Ts*.5 Ts*.5]; % The inital guess of time intervals
% u = [Vg  1 1 Io]';


ON = 1;
OFF = 0;

 Binary_for_DAB = [
ON OFF 
OFF OFF 
     ];

SW_OFF = ones(1,2).*10000000;

SW_ON = [M1_R,M2_R];

SW = [SW_OFF;SW_ON;SW_ON];




%}

%% The test case assignment
%%{

Numerical_Components = obj.Converter.Element_Properties;

for i = 1:1:size(Numerical_Components,1)
    eval([Numerical_Components{i,1} '=' 'Numerical_Components{i,2};'])
end


Resistances = obj.Converter.Switch_Resistors ;
SW = obj.Converter.Switch_Resistor_Values; % [SW_OFF; SW_ON; SW_ON]


%}




out = {};

k = 1;

% Index fixing to get off = 1 diode = 2 and ON = 3


new_state(new_state==2)=3; % Turn logic 2 (ON) to a 3
new_state(new_state==1)=2; % Turn logic 1 (DON) to a 2
new_state(new_state==-1)=1; % Turn logic 0 (OFF) to a -1


The_Codex = obj.Converter.Topology.Parser.Codex;


for i = 1:1:size(The_Codex,2)
    Binary(i,:) = new_state(The_Codex(:,i)); % Assign correct row
end

for j = 1:1:size(Binary,1)
    out{j} = SW(Binary(j),j);
end

for i = 1:1:size(Resistances,1)
    eval([Resistances{i} '=' 'out{i};'])
end

% [M1_R, M2_R, M3_R, M4_R, M5_R, M6_R, M7_R, M8_R] = deal(out{:});
% [M1_R, M2_R] = deal(out{:});
% [M1_R, M2_R, M3_R, M4_R, M5_R, M6_R, M7_R, M8_R,M9_R,M10_R, M11_R, M12_R] = deal(out{:});
% [M1_R, M2_R] = deal(out{:});

if isempty(obj.Converter.Topology.Parser.Asym)
    
    HtempAB(:,:,k) = eval(obj.Converter.Topology.Parser.HtempAB(:,:,1));
    HtempCD(:,:,k) = eval(obj.Converter.Topology.Parser.HtempCD(:,:,1));
    dependsAB(:,:,k) = eval(obj.Converter.Topology.Parser.dependsAB(:,:,1));
    savedCD(:,:,k) = eval(obj.Converter.Topology.Parser.savedCD(:,:,1));
    for j = 1:1:size(obj.Converter.Topology.Parser.DependentNames(:,1),1)
        DependentNames(j,k) = eval(obj.Converter.Topology.Parser.DependentNames{j,1});
    end
    for j = 1:1:size(obj.Converter.Topology.Parser.OutputNames(:,1),1)
        OutputNames(j,k) = eval(obj.Converter.Topology.Parser.OutputNames{j,1});
    end
    
    for k = 1:1:size(HtempAB,3)
        [A,B,C,D] = obj.Converter.Topology.Parser.loopfixAB_large(HtempAB(:,:,k),dependsAB(:,:,k),OutputNames(:,k),DependentNames(:,k));
        [C,D] = obj.Converter.Topology.Parser.loopfixCD_large(B,C,D,HtempCD(:,:,k),savedCD(:,:,k),DependentNames(:,k),obj.Converter.Topology.Parser.SortedTree(:,:,1),obj.Converter.Topology.Parser.SortedCoTree(:,:,1));
    end
    eigA = eig(A);
    
else
    
    A = eval(obj.Converter.Topology.Parser.Asym(:,:,1));
    B = eval(obj.Converter.Topology.Parser.Bsym(:,:,1));
    C = eval(obj.Converter.Topology.Parser.Csym(:,:,1));
    D = eval(obj.Converter.Topology.Parser.Dsym(:,:,1));
    eigA = eig(A);
    
end


% Sets new B and D matricies corresponding to the updated diode
% conduction indecies in the variable binary
B(:,contains(obj.Converter.Topology.Parser.ConstantNames,'VF'))=B(:,contains(obj.Converter.Topology.Parser.ConstantNames,'VF')).*(Binary==2)';
D(:,contains(obj.Converter.Topology.Parser.ConstantNames,'VF'))=D(:,contains(obj.Converter.Topology.Parser.ConstantNames,'VF')).*(Binary==2)';



end

