function [pass] = Run_Test_Case(selection,print)

%% Test Cases
% This code runs a test case of AURA and compares the RMS value of the
% result to a goal end value and shows the times that were used to
% complete the sequence

% Might eventually make a function to determine the test case wanted
% or to run all of them to run a final check






tic
switch selection
    
    case 'Test Primary Flyback'
    
      Assignment = tic;
        filename = 'EVAL_FLY_reduced.net';
        parse = NetListParse();
        
        Voltage = {'V1'
            'R1'
            'M1'
            'L2'};
        Current = {'V1'
            'R1'
            'M1'};
        
        parse.initialize(filename,Voltage,Current);
        parse.ABCD();
        
        top = SMPStopology();
        top.Parser = parse;
        
        conv = SMPSconverter();
        conv.Topology = top;
        
        simulator = SMPSim();
        
        Vg = 5;
        
        fs = 200e3;
        fs = 2.0367e+05;
        Ts = 1/fs;
        dt = 0.02/fs;
        dt = 0.1*Ts;
        
        ShortedL = 0.47382e-6;
        ShortedR = 326.9e-3;
        
        L4 = 1.0393e-6-ShortedL;  % Resonate indcutor
        L3 = 16.228e-6-ShortedL; % Magnetizing inductance
        
        R10 = 1000000; % Parallel to L4
        
        L1 = 1e-3; % Primary turns
        L2 = 1e-3; % Secondary turns
        C2 = 220e-9; % Resonate Cap
        C1 = 158e-6; % Output Cap
        
        
        
        R1 = 17;
        R2 = 390;
        R3 = 0.63521-ShortedR; % Inductor series resistance
        R9 = 0.001; % input resistance 
        
        M1_R = 0.054;
        M2_R = 0.0305;
        M3_R = 0.2834;
        
        M1_C = 166e-12;
        M2_C = 272e-12;
        M3_C = 28.47e-12;
        
        M3_C_Resist = 10000;
        M2_C_Resist = 10000;
        M1_C_Resist = 10000;
       
        
        u = [Vg 0.2 0.2]';
        Order  = [1 2];
        ts = [2.575e-6 Ts-2.575e-6]; % The inital guess of time intervals
        
        Numerical_Components = {'C1' C1
            'C2' C2
            'L1' L1
            'L2' L2
            'L3' L3
            'L3' L3
            'L4' L4
            'M1_C' M1_C
            'M2_C' M2_C
            'M3_C' M3_C
            'R1' R1
            'R2' R2
            'R3' R3
            'R9' R9
            'R10' R10
            'M3_C_Resist' M3_C_Resist
            'M2_C_Resist' M2_C_Resist
            'M1_C_Resist' M1_C_Resist
            };
        
        Switch_Resistors = {'M1_R'
            'M3_R'};
        
        ON = 1;
        OFF = 0;
        
        Binary_for_DAB = [
            ON OFF 
            OFF OFF
            ];
        
        SW_OFF = ones(1,2).*10000000;
        
        SW_ON = [M1_R,M3_R];
        
        SW = [SW_OFF;SW_ON;SW_ON];
    
        
        
        
    
    case 'EVAL-CN0342-EB1Z'
        
        Assignment = tic;
        filename = 'EVAL_FLY.net';
        parse = NetListParse();
        
        Voltage = {'V1'
            'R1'
            'M1'
            'L2'};
        Current = {'V1'
            'R1'
            'M1'};
        
        parse.initialize(filename,Voltage,Current);
        parse.ABCD();
        
        top = SMPStopology();
        top.Parser = parse;
        
        conv = SMPSconverter();
        conv.Topology = top;
        
        simulator = SMPSim();
        
        Vg = 5;
        
        fs = 200e3;
        fs = 2.0367e+05;
        Ts = 1/fs;
        dt = 0.02/fs;
        dt = 0.1*Ts;
        
        ShortedL = 0.47382e-6;
        ShortedR = 326.9e-3;
        
        L4 = 1.0393e-6-ShortedL;  % Resonate indcutor
        L3 = 16.228e-6-ShortedL; % Magnetizing inductance
        
        R10 = 1000000; % Parallel to L4
        
        L1 = 1e-3; % Primary turns
        L2 = 1e-3; % Secondary turns
        C2 = 220e-9; % Resonate Cap
        C1 = 158e-6; % Output Cap
        
        
        
        R1 = 17;
        R2 = 390;
        R3 = 0.63521-ShortedR; % Inductor series resistance
        R9 = 0.001; % input resistance 
        
        M1_R = 0.054;
        M2_R = 0.0305;
        M3_R = 0.2834;
        
        M1_C = 166e-12;
        M2_C = 272e-12;
        M3_C = 28.47e-12;
        
        M3_C_Resist = 10000;
        M2_C_Resist = 10000;
        M1_C_Resist = 10000;
       
        
        u = [Vg 0.2 0.2 0.2]';
        Order  = [1 2];
        ts = [2.575e-6 Ts-2.575e-6]; % The inital guess of time intervals
        
        Numerical_Components = {'C1' C1
            'C2' C2
            'L1' L1
            'L2' L2
            'L3' L3
            'L3' L3
            'L4' L4
            'M1_C' M1_C
            'M2_C' M2_C
            'M3_C' M3_C
            'R1' R1
            'R2' R2
            'R3' R3
            'R9' R9
            'R10' R10
            'M3_C_Resist' M3_C_Resist
            'M2_C_Resist' M2_C_Resist
            'M1_C_Resist' M1_C_Resist
            };
        
        Switch_Resistors = {'M1_R'
            'M2_R'
            'M3_R'};
        
        ON = 1;
        OFF = 0;
        
        Binary_for_DAB = [
            ON OFF OFF 
            OFF OFF OFF
            ];
        
        SW_OFF = ones(1,3).*10000000;
        
        SW_ON = [M1_R,M2_R,M3_R];
        
        SW = [SW_OFF;SW_ON;SW_ON];
        
        
    case 'HDSC'
        
        
        Assignment = tic;
        filename = 'HDSC_AURA.net';
        parse = NetListParse();
        
        Voltage = {'V1'
            'R1'
            'M1'};
        Current = {'V1'
            'R1'
            'M1'};
        
        parse.initialize(filename,Voltage,Current);
        parse.ABCD();
        
        top = SMPStopology();
        top.Parser = parse;
        
        conv = SMPSconverter();
        conv.Topology = top;
        
        simulator = SMPSim();
        
        
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
        
        
        
        Numerical_Components = {'V1' Vg
            'C1' C1
            'C2' C2
            'C3' C3
            'C4' C4
            'L1' L1
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
            'R4' R4
            'R5' R5};
        
        Switch_Resistors = {'M1_R'
            'M2_R'
            'M3_R'
            'M4_R'
            'M5_R'
            'M6_R'
            'M7_R'
            'M8_R'};
        

    
    case 'BuckBoost'
        
        Assignment = tic;
        filename = 'Buck_Boost.net';
        parse = NetListParse();
        parse.initialize(filename);
        parse.ABCD();
        
        top = SMPStopology();
        top.Parser = parse; 
        
        conv = SMPSconverter();
        conv.Topology = top;
        
        simulator = SMPSim();
        
        
        
        
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
         
         
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'DAB'
        Assignment = tic;
        filename = 'DAB_Resistors_Cap.net';
        parse = NetListParse();
        parse.initialize(filename);
        parse.ABCD();
        
        top = SMPStopology();
        top.Parser = parse;
        
        conv = SMPSconverter();
        conv.Topology = top;
        
        simulator = SMPSim();
        
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
        
        ts = [primary_dead PS secondary_dead Power primary_dead PS secondary_dead Power ];
        
        % This is the the input voltage and the all of the body diode forward
        % voltages
        u = [Vg 1 1 1 1 1 1 1 1]';
        
        TestparseWaveform = false;
        
        ON = 1;
        OFF = 0;
        
        
        SW_OFF = ones(1,8).*10000000;
        
        SW_ON = [0.05 0.05 0.05 0.05 0.0015 0.0015 0.0015 0.0015];
        
        SW = [SW_OFF;SW_ON;SW_ON];
        
        Binary_for_DAB = [
            
        OFF OFF OFF OFF ON OFF OFF ON %primary sw
        OFF ON ON OFF ON OFF OFF ON % phase shift
        OFF ON ON OFF OFF OFF OFF OFF % secondary sw
        OFF ON ON OFF OFF ON ON OFF % Reverse power
        OFF OFF OFF OFF OFF ON ON OFF % primary sw
        ON OFF OFF ON OFF ON ON OFF % phase shift
        ON OFF OFF ON OFF OFF OFF OFF % secondary sw
        ON OFF OFF ON ON OFF OFF ON]; % POWER
    
    
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
        'R3' R3 };
    
    
    Switch_Resistors = {'M1_R'
        'M2_R'
        'M3_R'
        'M4_R'
        'M5_R'
        'M6_R'
        'M7_R'
        'M8_R'};
    
    
    THE_KEY = [
        33.8525
        33.8525
        0.9507
        0.9507
        0.9507
        0.0021
        0.2748
        33.8525
        33.8525
        0.9507
        1.3229];
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'Ind_Dickson'
        
        Assignment = tic;
        filename = 'Induct_Load_Dickson.net';
        parse = NetListParse();
        parse.initialize(filename);
        parse.ABCD();
        
        top = SMPStopology();
        top.Parser = parse;
        
        conv = SMPSconverter();
        conv.Topology = top;
        
        simulator = SMPSim();
        
        
        Vg = 48;
        Iload = 5;
        fs = 1e6;
        Ts = 1/fs;
        dt = Ts/100;
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
        
        
        u = [Vg  1 1 1 1 1 1 1 1 1 1 1 1 Iload]';
        Order  = [1 2 3 4];
        ts = [Ts*.5-dt dt Ts*.5-dt dt]; % The inital guess of time intervals
        
        
        
        ON = 1;
        OFF = 0;
        
        
        % Binary_for_DAB = [
        % OFF ON OFF ON OFF ON OFF ON ON OFF OFF ON
        % ON OFF ON OFF ON OFF ON OFF OFF ON ON OFF
        %      ];
        
        Binary_for_DAB = [
            OFF ON OFF ON OFF ON OFF ON ON OFF OFF ON
            OFF OFF OFF OFF OFF OFF OFF OFF OFF OFF OFF OFF
            ON OFF ON OFF ON OFF ON OFF OFF ON ON OFF
            OFF OFF OFF OFF OFF OFF OFF OFF OFF OFF OFF OFF
            ];
        
        SW_OFF = ones(1,12).*10000000;
        
        SW_ON = [M1_R,M2_R,M3_R,M4_R,M5_R,M6_R,M7_R,M8_R,M9_R,M10_R,M11_R,M12_R];
        
        SW = [SW_OFF;SW_ON];
        
        Numerical_Components = {'V1' Vg
            'C1' C1
            'C2' C2
            'C3' C3
            'C4' C4
            'C5' C5
            'C6' C6
            'C7' C7
            'C8' C8
            'L1' L1
            'M1_C' M1_C
            'M2_C' M2_C
            'M3_C' M3_C
            'M4_C' M4_C
            'M5_C' M5_C
            'M6_C' M6_C
            'M7_C' M7_C
            'M8_C' M8_C
            'M9_C' M9_C
            'M10_C' M10_C
            'M11_C' M11_C
            'M12_C' M12_C
            'R1' R1
            'R2' R2
            'R3' R3
            'R4' R4
            'R5' R5
            'R6' R6
            'R7' R7
            'R8' R8};
        
        Switch_Resistors = {'M1_R'
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
            'M12_R'};
        
        THE_KEY = [
            4.2401
            8.4303
            8.5074
            8.5074
            42.0640
            36.0787
            30.0397
            24.0001
            8.5074
            8.5074
            8.4303
            4.2404
            17.9611
            11.9216
            5.9377
            5.8051
            4.1757
            4.1795
            5.0001
            4.1760
            4.1791];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'MR_Buck'
        
        Assignment = tic;
        filename = 'MR_Buck.net';
        parse = NetListParse();
        parse.initialize(filename);
        parse.ABCD();
        
        top = SMPStopology();
        top.Parser = parse;
        
        conv = SMPSconverter();
        conv.Topology = top;
        
        simulator = SMPSim();
        
        Vg = 12;
        L1 = 1e-6; %L
        C1 = 1*10^-6; %Cout
        L2 = 0.01e-9; % Resonate inductor
        fs = 1e6;
        Ts = 1/fs;
        V = 3;
        Io = 0.1; % was 1
        M1_C = 1.5e-9; % CHS
        M2_C = 1.5e-9; % LHS

        M1_R_ON = 0.002; % ronHS
        M2_R_ON = 0.002; % ronLS

        Order = [1 2]; % The order that the states must go in after being parsed
        
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
        
        SW = [SW_OFF; SW_ON; SW_ON];
        
          Numerical_Components = {
            'C1' C1
            'L1' L1
            'L2' L2
            'M1_C' M1_C
            'M2_C' M2_C
            };
        
        
        Binary_for_DAB = [
            ON OFF
            OFF OFF];
        
        
        Switch_Resistors = {'M1_R'
            'M2_R'};
        
        THE_KEY = [
        11.1089
        2.4463
        6.2030
        0.7540
        7.9659];
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    case 'AC Fly'
        
        Assignment = tic;
        filename = 'AC Flyback.net';
        parse = NetListParse();
        parse.initialize(filename);
        parse.ABCD();
        
        top = SMPStopology();
        top.Parser = parse;
        
        conv = SMPSconverter();
        conv.Topology = top;
        
        simulator = SMPSim();
        
        Vg = 5;
        
        fs = 1e6;
        Ts = 1/fs;
        dt = 0.02/fs;
        L1 = 1e-3; % Primary turns
        L2 = 16e-3; % Secondary turns
        L3 = 0.98648e-6;  % Magnetizing inductance
        r_on = 16e-3;
        Cds = 150e-12;
        L4 = 50e-9; % Resonate indcutor
        C2 = 111.91e-9; % Resonate Cap
        C1 = 10e-6; % Output Cap
        R1 = 133.3333;
        R2 = 0.001;
        R4 = 0.001;
        R3 = 1000;
        I1 = 0.15;
        [M1_R,M2_R,M3_R] = deal(r_on);
        [M1_C,M2_C,M3_C] = deal(Cds);
        M3_C_Resist = 10000;
        M2_C_Resist = 10000;
        M1_C_Resist = 10000;
       
        
        u = [Vg 1 1 1]';
        Order  = [1 2 3 4];
        ts = [Ts*.53-dt dt Ts*.47 dt]; % The inital guess of time intervals
        
        Numerical_Components = {'C1' C1
            'C2' C2
            'L1' L1
            'L2' L2
            'L3' L3
            'L4' L4
            'M1_C' M1_C
            'M2_C' M2_C
            'M3_C' M3_C
            'R1' R1
            'R2' R2
            'R3' R3
            'R4' R4
            'I1' I1
            'M3_C_Resist' M3_C_Resist
            'M2_C_Resist' M2_C_Resist
            'M1_C_Resist' M1_C_Resist
            };
        
        Switch_Resistors = {'M1_R'
            'M2_R'
            'M3_R'};
        
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
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    case 'COMPEL Buck'
        
        Assignment = tic;
        filename = 'Buck_Qual.net';
        parse = NetListParse();
        parse.initialize(filename);
        parse.ABCD();
        
        top = SMPStopology();
        top.Parser = parse;
        
        conv = SMPSconverter();
        conv.Topology = top;
        
        simulator = SMPSim();
        %{
            Ls = 1e-6;
    Ron = 0.0015;
    Vg = 12;
    Vout = 5;
    % R_L should be around 0.16 for experimetnal stuff
    I_out = 3;
    Cout = (4*100+680)*10^-6;
    R_out = 1;
    R_L = 0.16;
    Coss = 1620e-12;
    Vfwd = 1;
    
    
    D = 5/12;
    fs_adj = 200e3;
    Ts = 1/fs_adj;
    dead = 0.001;
    M2_D = (1-D)-dead;
    M1_D = D-dead;
    M2_delay = Ts*(D+dead);
    M1_delay = 0;
    %}
        
        
            D = 5/12;
    fs_adj = 200e3;
    Ts = 1/fs_adj;
    dead = 0.002;
        
        
        
        Vg = 12;
        
        fs = 0.2e6;
        Ts = 1/fs;
        dt = 0.02/fs;
        L1 = 6e-6; 
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
       
        
        u = [Vg 1 1]';
     %   Order  = [1 2 3 4];
        ts = [Ts*(D-dead) Ts*(dead) Ts*(1-D-dead) Ts*(dead)]; % The inital guess of time intervals
        
        Numerical_Components = {'C1' C1
            'L1' L1
            'M1_C' M1_C
            'M2_C' M2_C
            'R1' R1
            'R2' R2
            'R3' R3
            'M2_C_Resist' M2_C_Resist
            'M1_C_Resist' M1_C_Resist
            };
        
        Switch_Resistors = {'M1_R'
            'M2_R'};
        
        ON = 1;
        OFF = 0;
        
        Binary_for_DAB = [
            ON OFF 
            OFF OFF 
            OFF ON 
            OFF OFF 
            
            ];
        
        SW_OFF = ones(1,2).*10000000;
        
        SW_ON = [M1_R,M2_R];
        
        SW = [SW_OFF;SW_ON;SW_ON];
        
        
        
    otherwise
        fprintf('Pick a test case that actually exists please\n Number of times this has happend: %d \n',ceil(rand*1000)) % That's embarrassing
        return
end

conv.ts = ts;
conv.u = u;
conv.order = Order;
conv.Element_Properties = Numerical_Components;
conv.Switch_Resistors = Switch_Resistors;
conv.Switch_Resistor_Values = SW;




%{
% Switch Resistors Contains the same order of resistors as was given
% in the Binary matrix provided by the user to determine when states
% are on or off
for i = 1:1:size(Switch_Resistors,1)
    eval([Switch_Resistors{i} '=' 'out{i}']);
end
%}



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
        for i = 1:1:size(Switch_Resistors,1)
            eval([Switch_Resistors{i} '=' 'out{i};'])
        end
        
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
        [parse.eigA(:,k)] = eig(parse.Anum(:,:,k));
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
        
        for i = 1:1:size(Switch_Resistors,1)
            eval([Switch_Resistors{i} '=' 'out{i};']);
        end
        
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
Xss_Aug=simulator.SS_Soln_Aug();

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


parse.find_diode_new(Order,Binary_for_DAB);
iterations = 100;
Optimization = SMPSOptim;
Optimization.Simulator = simulator;
Optimization.opptimization_loop;
% simulator.Three_tier_diode_correct(iterations,print)

[~, ~, y] = simulator.SS_WF_Reconstruct();

StateNumbers = parse.StateNumbers;

State_RMS = rms(y(StateNumbers,:),2);




difference = State_RMS-THE_KEY;

if print
    fprintf('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n')
    fprintf('Test_Case %s Results\n',selection)
    fprintf('------------------------------\n')
    fprintf('%15s %12s \n','State Variable','Difference')
    for i  = 1:length(difference)
        fprintf('%15s %12.8f \n',top.stateLabels{i},difference(i))
    end
    fprintf('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n')
end

if sum(difference>0.5)>0
    pass = false;
    fprintf('Better luck next time')
else
    pass = true;
end


return















