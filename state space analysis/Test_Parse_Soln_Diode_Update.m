% This script tests the combination of state space analysis and circuit
% parsing 

% Currently trying to get Buck converter to work 
% Then work on making it general for other tests

%     _   _   _  ____    _    
%    / \ | | | |/ _  |  / \   
%   / _ \| | | | (_| | / _ \  
%  / ___ | |_| |> _  |/ ___ \ 
% /_/   \_\___//_/ |_/_/   \_\


%% User Input

clear
close all
tic
%% Boost converter
%{
V = 5;
L1 = 230e-9; %L
C1 = 4040e-9; %Cout
fs = 2e6;
Ts = 1/fs;
Vg = 1.8;
Io = 10; % was 1
M1_C = 3.4874e-10; % CHS
M2_C = 3.4874e-10; % LHS
D1_C = 3.4874e-10; % LHS

R1 =  5; % Output resistor
dt = Ts/100;%5e-10;
Vdr = 5;
M1_R_ON = .05; % ronHS
M2_R_ON = .05; % ronLS
D1_R_D = .05; % ronLS
M1_R_D = .05; % ronHS
M2_R_D = .05; % ronLS
[D1_R_OFF,M1_R_OFF,M2_R_OFF] = deal(100000000); % ronLS
Order = [2 1 4 1]; % The order that the states must go in after being parsed
ts = [Ts*.5-dt dt Ts*.5-dt dt]; % The inital guess of time intervals
u = [Vg 1 1]';
%}

%% Buck Converter
%{
Vg = 5;
L1 = 230e-9; %L
C1 = 4040e-9; %Cout
fs = 2e6;
Ts = 1/fs;
V = 1.8;
Io = 1; % was 1
M1_C = 3.4874e-10; % CHS
M2_C = 3.4874e-10; % LHS
D1_C = 3.4874e-10; % LHS

R2 = 0.5;
R1 =  .01; % Rl also known as R1 (one or L (lowercase))
dt = Ts/100;%5e-10;
Vdr = 5;
M1_R_ON = .05; % ronHS
M2_R_ON = .05; % ronLS
D1_R_D = .05; % ronLS
M1_R_D = .05; % ronHS
M2_R_D = .05; % ronLS
[D1_R_OFF,M1_R_OFF,M2_R_OFF] = deal(100000000); % ronLS

%Order = [2 4]; % The order that the states must go in after being parsed
%ts = [Ts*.5 Ts*.5]; % The inital guess of time intervals


Order = [2 1 4 1]; % The order that the states must go in after being parsed
ts = [Ts*.5-dt dt Ts*.5-dt dt]; % The inital guess of time intervals
u = [Vg  1 1 Io]';
%}

%% Flyback Converter
%{
Vg = 12;
L1 = 1e-3;
L2 = 4e-3;
L3 = 0.5e-3;
C1 = 4040e-9; %Cout
fs = 2e6;
Ts = 1/fs;
V = 1.8;
Io = 10; % was 1
M1_C = 3.4874e-10; % CHS
M2_C = 3.4874e-10; % LHS
D1_C = 3.4874e-10; % LHS

R1 =  24; % Rl also known as R1 (one or L (lowercase))
dt = Ts/1000;%5e-10;
Vdr = 5;
M1_R_ON = .05; % ronHS
M2_R_ON = .05; % ronLS
D1_R_D = .05; % ronLS
M1_R_D = .05; % ronHS
M2_R_D = .05; % ronLS
[D1_R_OFF,M1_R_OFF,M2_R_OFF] = deal(100000000); % ronLS

Order = [2 1 4 1]; % The order that the states must go in after being parsed
ts = [Ts*.66-dt dt Ts*.34-dt dt]; % The inital guess of time intervals


%Order = [2 1 4 1]; % The order that the states must go in after being parsed
%ts = [Ts*.5-dt dt Ts*.5-dt dt]; % The inital guess of time intervals
u = [Vg  1 1]';
%}

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

SW = [SW_OFF;SW_ON];


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

%{
% The Normal ON
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

%%{
%The FLIPED one
Binary_for_DAB = [

    OFF OFF OFF OFF ON OFF OFF ON %primary sw
    OFF ON ON OFF ON OFF OFF ON % phase shift
    OFF ON ON OFF OFF OFF OFF OFF % secondary sw
    OFF ON ON OFF OFF ON ON OFF % Reverse power
    OFF OFF OFF OFF OFF ON ON OFF % primary sw
    ON OFF OFF ON OFF ON ON OFF % phase shift
    ON OFF OFF ON OFF OFF OFF OFF % secondary sw
    ON OFF OFF ON ON OFF OFF ON]; % POWER




for k = 1:1:size(Binary_for_DAB,1)
        
        out = {};
        
        for j = 1:1:size(Binary_for_DAB,2)
            out{j} = SW(Binary_for_DAB(k,j)+1,j);
        end
         [M1_R, M2_R, M3_R, M4_R, M5_R, M6_R, M7_R, M8_R] = deal(out{:});

end

M1_C_VF = 0;
M2_C_VF = 0;
M3_C_VF = 0;
M4_C_VF = 0;
M5_C_VF = 0;
M6_C_VF = 0;
M7_C_VF = 0;
M8_C_VF = 0;



Numerical_Components = {'V1' Vg
    'C1' C1
    'C2' C2
    'L1' L1
    'L2_T' L2
    'L3_T' L3
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
    'M1_C_VF' M1_C_VF
    'M2_C_VF' M2_C_VF
    'M3_C_VF' M3_C_VF
    'M4_C_VF' M4_C_VF
    'M5_C_VF' M5_C_VF
    'M6_C_VF' M6_C_VF
    'M7_C_VF' M7_C_VF
    'M8_C_VF' M8_C_VF
    'M1_R' M1_R
    'M2_R' M2_R
    'M3_R' M3_R
    'M4_R' M4_R
    'M5_R' M5_R
    'M6_R' M6_R
    'M7_R' M7_R
    'M8_R' M8_R };



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
%dead = 0.001/fs;
dead = 0.1/fs;
ts = [Ts*D-dead dead Ts*(1-D)-dead dead];

u = [Vg 1 1]';


TestparseWaveform = false;

ON = 1;
OFF = 0;


SW_OFF = ones(1,2).*10000000;

SW_ON = [M1_R_ON M2_R_ON];

SW = [SW_OFF;SW_ON];


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
%}

%% MR-Buck Converter
%%{
Vg = 12;
L1 = 1e-6; %L
C1 = 1*10^-6; %Cout
L2 = 0.01e-9; % Resonate inductor
L3 = 0.01e-9; % Resonate inductor
fs = 1e6;
Ts = 1/fs;
V = 3;
Io = 0.1; % was 1
M1_C = 1.5e-9; % CHS
M2_C = 1.5e-9; % LHS
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


Order = [1 2]; % The order that the states must go in after being parsed
% ts = [Ts*.5-dt dt Ts*.5-dt dt]; % The inital guess of time intervals
% u = [Vg  1 1 Io]';


% Order = [1 2 3 4 5 6 7 8]; % The order that the states must go in after being parsed
% ts = [Ts*.2376 Ts*.0277 Ts*.23 Ts*.00023 Ts*.2376 Ts*.02775 Ts*.23 Ts*.00023];



%[D1_R,M1_R,M2_R,M3_R,M4_R,M5_R,M6_R,M7_R,M8_R] = deal(100000000); % ronLS

D = 5/12;
%dead = 0.001/fs;
dead = 0.01/fs;
% ts = [25e-9 25e-9 25e-9 25e-9];
ts = [Ts/4 3*Ts/4];
u = [Vg 1 1 Io]';


TestparseWaveform = false;

ON = 1;
OFF = 0;


SW_OFF = ones(1,2).*10000000;

SW_ON = [M1_R_ON M2_R_ON];

SW = [SW_OFF;SW_ON];


%Out = [M1_R,M2_R,M3_R,M4_R,M5_R,M6_R,M7_R,M8_R];



% Binary_for_DAB = [
%     ON ON
%     ON OFF
%     ON ON
%     OFF ON];


Binary_for_DAB = [
    ON OFF
    OFF OFF];




%}




% Select .net file
% filename = 'DAB_Resistors_Cap.net';
% filename = 'HDSC_AURA.net';
% filename = 'Buck_Qual.net';
%   filename = 'MR_Buck.net';
filename = 'Buck_io.net';
% Current options for filename:
% Boost.net
% Buck.net
% Buck2.net This is a buck with a series inductor resistance and current source output 
% Dickson.net
% Flyback.net
% Forward.net
% DAB.net
% DAB_no_RL.net
% Buck_Qual.net Synchronous buck with series inductor resistance

% Select test case file
testcase = 'TEST_ABCD_Buck_Diode';
% Current options for testcase:
% TEST_ABCD_Boost
% TEST_ABCD_Buck
% TEST_ABCD_Buck2
% TEST_ABCD_Flyback
% TEST_ABCD_Forward
% TEST_ABCD_Dickson

% Set Voltage and Current Nodes to add
Voltage = {'V1'
    'R1'
    'M1'};
Current = {'V1'
    'R1'
    'M1'};
% Change Voltage and Current based on desired output measurements (C and D
% matricies). Voltage and Current should be of type Cell 
% Example:
% Voltage = {'V1'
%     'M1'
%     'L3'};
% 
% Current = {'C1'
%     'D2'
%     'R3'};

%%% Note All Swithch, Inductor, and Capaictor elements will be included in
%%% the C and D matrix whether they are placed in the above cell array or
%%% not.



%% Run functions
tic
parse = NetListParse();
%parse.initialize(filename,Voltage,Current,Numerical_Components);
parse.initialize(filename,Voltage,Current);
parse.ABCD();
%  load('D:\GitHub\AURA\DAB_PARSE.mat')


% if TestparseWaveform
% testfun = str2func(testcase);
% testfun(parse);
% end

%% DC code (set up converter and topology classes)
toc
top = SMPStopology();
top.Parser = parse;

conv = SMPSconverter();
conv.Topology = top;
conv.ts = ts;
conv.u = u;
conv.order = Order;
%  parse.find_diode(Order);

simulator = SMPSim();

%% To update and test
    % taken from TEST function

A = parse.Asym;
B = parse.Bsym;
C = parse.Csym;
D = parse.Dsym;

SortedTree = parse.SortedTree;
SortedCoTree = parse.SortedCoTree;

if isempty(A)
    for k = 1:1:size(Binary_for_DAB,1)
        
        out = {};
        
        for j = 1:1:size(Binary_for_DAB,2)
            out{j} = SW(Binary_for_DAB(k,j)+1,j);
        end
        % [M1_R, M2_R, M3_R, M4_R, M5_R, M6_R, M7_R, M8_R] = deal(out{:});
         [M1_R, M2_R] = deal(out{:});
        
        
        HtempAB(:,:,k) = eval(parse.HtempAB(:,:,1));
        HtempCD(:,:,k) = eval(parse.HtempCD(:,:,1));
        dependsAB(:,:,k) = eval(parse.dependsAB(:,:,1));
        savedCD(:,:,k) = eval(parse.savedCD(:,:,1));
        for j = 1:1:size(parse.DependentNames(:,1),1)
            DependentNames(j,k) = eval(parse.DependentNames{j,1});
        end
        for j = 1:1:size(parse.OutputNames(:,1),1)
            OutputNames(j,k) = eval(parse.OutputNames{j,1});
        end
    end
    for k = 1:1:size(HtempAB,3)
        [A,B,C,D] = parse.loopfixAB_large(HtempAB(:,:,k),dependsAB(:,:,k),OutputNames(:,k),DependentNames(:,k));
        [C,D] = parse.loopfixCD_large(B,C,D,HtempCD(:,:,k),savedCD(:,:,k),DependentNames(:,k),SortedTree(:,:,1),SortedCoTree(:,:,1));
        parse.Anum(:,:,k)=A;
        parse.Bnum(:,:,k)=B;
        parse.Cnum(:,:,k)=C;
        parse.Dnum(:,:,k)=D;
    end
    
    % Set all diode forward voltages to be off
    B = parse.Bnum;
    B(:,contains(parse.ConstantNames,'VF'),:)=0;
    parse.Bnum = B;
    D = parse.Dnum;
    D(:,contains(parse.ConstantNames,'VF'),:)=0;
    parse.Dnum = D;
    
elseif ~isempty(parse.Anum)
    fprintf('Confirm Plecs used');
else
    for k = 1:1:size(Binary_for_DAB,1)
        
        out = {};
        
        for j = 1:1:size(Binary_for_DAB,2)
            out{j} = SW(Binary_for_DAB(k,j)+1,j);
        end
        
        %[M1_R, M2_R, M3_R, M4_R, M5_R, M6_R, M7_R, M8_R] = deal(out{:});
        [M1_R, M2_R] = deal(out{:});
        
        parse.Anum(:,:,k) = eval(A(:,:,1));
        parse.Bnum(:,:,k) = eval(B(:,:,1));
        parse.Cnum(:,:,k) = eval(C(:,:,1));
        parse.Dnum(:,:,k) = eval(D(:,:,1));
        [parse.eigA(:,k)] = eig(parse.Anum(:,:,k));
    end
     % Set all diode forward voltages to be off
    B = parse.Bnum;
    B(:,contains(parse.ConstantNames,'VF'),:)=0;
    parse.Bnum = B;
    D = parse.Dnum;
    D(:,contains(parse.ConstantNames,'VF'),:)=0;
    parse.Dnum = D;
end

simulator.loadTestConverter2(conv);
simulator.eigA = parse.eigA;
simulator.binary = Binary_for_DAB;
Xss = simulator.SS_Soln();

%% Reconstruction of Dependent variables
parse.StateVarIndex();
simulator.CorrectXs();


for i = 1:1:length(parse.StateNumbers)
    if strcmp(parse.OutputNamesCD{parse.StateNumbers(i),1}(1,end),'A')
        top.stateLabels(end+1,1) = strcat('I_{', parse.StateNames(i,1),'} (A)');
        top.stateLabels_Opp(end+1,1) = strcat('V_{', parse.StateNames(i,1),'} (V)');
    else
        top.stateLabels(end+1,1) = strcat(parse.OutputNamesCD{parse.StateNumbers(i),1}(1,end),'_{', parse.StateNames(i,1),'} (V)');
        top.stateLabels_Opp(end+1,1) = strcat('I_{', parse.StateNames(i,1),'} (A)');
    end
end
    
    
    % Dont need anymore due to fix of SS_Soln.m
% Dependent variables are calculated with independent variables
%{
Xss(end+1,:) = zeros(size(parse.DependentNames,1),size(Xss,2));


% For now diode is 8th row in C and D will char match dependent variables
% to find where they are in output and calculate
for i = 2:1:size(Xss,2)
    Xss(size(OutputNames,1)+1:end,i) = simulator.Cw(8,:,i-1)*Xss(:,i)+simulator.Dw(8,:,i-1)*u;
    Xss(:,1) = Xss(:,end);
end
%}

%% Attempt to solve for t directly

start_value = simulator.Xs(:,2);
end_value = simulator.Xs(:,2);
end_value(1) = 5;
end_value(4) = 0;
time_pos = 2;
% t = simulator.deadtimecalc(start_value,end_value,time_pos);

%% adjustDiodeCond


% User input to define output voltage/Cap
% identify output cap/component

optimize = SMPSOptim;
optimize.Simulator = simulator;
Vopos = 2;
Xs = Xss;

optimize.dead_time_intervals = [2,4];
optimize.dead_time_states = [4,2];
optimize.dead_time_goals = [0,0];

optimize.Vo_index = 10 ;
optimize.Vo_ideal_value = V;
optimize.Perturb1_index = [1 5]; % Used for state sensitivity power time
optimize.Perturb2_index = [3 7]; % Used for state sensitivity power time
iterations = 10;
parse.find_diode_new(Order,Binary_for_DAB);
% check = simulator.VfwdIrev();

% simulator.DAB_Mesh;
optimize.opptimization_loop;

%check = simulator.VfwdIrev();

% Have two conditions to check physical validity of diodes
% Check to see if make sure they have no negative current (should be blocking)
% Check to see if the voltage exceeds the forward voltage of the diode.
% These do not need to know what the state of anything is wheter fet or
% power diode

return
%{
simulator.Baxter_adjustDiodeCond();




Voerr = mean(Xs(Vopos,:)) - V; % Error on output
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
%}


