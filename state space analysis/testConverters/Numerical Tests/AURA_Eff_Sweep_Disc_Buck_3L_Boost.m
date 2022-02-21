function [stick] = AURA_Eff_Sweep_Disc_Buck_3L_Boost(X)

%     _   _   _  ____    _
%    / \ | | | |/ _  |  / \
%   / _ \| | | | (_| | / _ \
%  / ___ | |_| |> _  |/ ___ \
% /_/   \_\___//_/ |_/_/   \_\



% Current konwn issues use FETs for diodes. (Just have them never turn
% on)

try
    
    % Place in the filename
    filename = 'Buck_Boost_Vout.net'; % Place netlist filename here that you want to run
    
    % Can also place component values here to get their voltage and
    % current output waveforms such as:
    
    Voltage = {'V1'
        'V2'
        'V3'
        'R7'
        };
    Current = {'V1'
        'V2'
        'V3'
        'R7'
        };
    
    
    %% Determining Coss and ron based on w
    
    
    
    
    penalty  = 0;
    adjust = [ones(1,6), 1, 1e-6, 100, 100] ;
    X = X./adjust;
    
    FET_selection = X(1:6);
    Select_L = X(end-3);
    fs = X(end-2);
    dt1 = X(end-1);
    dt2 = X(end);
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
    %{
    switch Select_L
        
        case 1
            %XEL3515-151
            [R1_L, R2_L, C_L, k1, k2, k3, k4, k5] = deal(57, 4.80E-03, 1.7, 1.87E-05, 3.41E-02, 0.15, 4.81E-04, 2.05E-7);
            L_area = 3.2*3.5;
            
        case 2
            %XEL3520-131
            [R1_L, R2_L, C_L, k1, k2, k3, k4, k5] = deal(99, 3.50E-03, 2.0, 4.73E-05, 5.50E-02, 0.133, 1.32E-03, 4.50E-06);
            L_area = 3.2*3.5;
            
        case 3
            %XGL4018-121
            [R1_L, R2_L, C_L, k1, k2, k3, k4, k5] = deal(2, 2.00E-03, 3.4, 1.0E-06, 0.017, 0.12, 1.0E-06, 1.0E-06);
            L_area = 4*4;
            
        case 4
            %XAL7020-151
            [R1_L, R2_L, C_L, k1, k2, k3, k4, k5] = deal(10, 0.0019, 6.9, 2.00E-05, 0.020,  0.15,  1.00E-04, 5.00E-06);
            L_area = 7.5*7.5;
            
            
        otherwise
            
            stick = 100;
            return
            
    end
    %}

     List_of_Inductors = [
        57        0.0048       1.7      1.87e-05     0.0341       0.15       0.000481    0.000205   12.2275
        99        0.0035         2      4.73e-05      0.055      0.133        0.00132     4.5e-06   12.2275
        19         0.005       4.6         5e-06      0.017       0.12          1e-06       3e-06   18.49
        2         0.002       3.4         1e-06      0.017       0.12          1e-06       1e-06    18.49
        11        0.0015      3.07         1e-06     0.0223       0.13          1e-06       1e-06   18.49
        6        0.0018       4.5         1e-06      0.022       0.15          1e-06       1e-06    18.49
        5       0.00215       4.3         5e-05      0.016       0.16          0.001       5e-06    31.1264
        66.4       0.00153      5.62         1e-06      0.184       0.14          0.001        0.01 31.1264
        33.7       0.00148       5.7         1e-06      0.046       0.14        0.00058        0.01 31.1264
        42          0.17       100       1.8e-06      0.016        0.1         0.0005     9.8e-06   52.82
        42          0.17       120      1.88e-06     0.0297      0.122        0.00119     9.8e-06   52.82
        39.6          0.17       145      3.28e-06      0.024      0.152        0.00151     9.8e-06 52.82
        10        0.0019       6.9         2e-05       0.02       0.15         0.0001       5e-06   64
        20       0.00115      7.02         4e-05      0.022       0.16         0.0005       6e-06   64
        19       0.00075      5.31         1e-06       0.02       0.16         0.0199       1e-06   61.6
        472         0.001       835      1.73e-05      0.115       0.18        0.00285    9.38e-06  54.229
        52       0.00039        90      1.32e-05      0.012      0.122       -0.00101    3.39e-06   71.4
        48       0.00039        75      1.61e-05      0.018      0.152       -0.00181    3.39e-06   71.4
        50       0.00039       115      4.52e-06      0.013      0.222       -0.00165    3.39e-06   71.4
        52       0.00047        90      1.32e-05      0.012      0.122       -0.00101    3.39e-06   71.4
        48       0.00047        75      1.61e-05      0.018      0.152       -0.00181    3.39e-06   71.4
        50       0.00047       115      4.52e-06      0.013      0.222       -0.00165    3.39e-06   71.4
        52       0.00055        90      1.32e-05      0.012      0.122       -0.00101    3.39e-06   71.4
        48       0.00055        75      1.61e-05      0.018      0.152       -0.00181    3.39e-06   71.4
        50       0.00055       115      4.52e-06      0.013      0.222       -0.00165    3.39e-06   71.4
        39       0.00048        86      1.32e-05      0.013      0.121         -0.001    3.39e-06   83.2
        41       0.00048        90      1.32e-05      0.014      0.141         -0.001    3.39e-06   83.2
        35       0.00048       110      1.32e-05      0.014      0.171         -0.001    3.39e-06   83.2
        43       0.00029        78      1.32e-05      0.013      0.121       -0.00101    3.39e-06   83.2
        41         0.029       103      1.32e-05      0.015      0.151       -0.00101    3.39e-06   83.2
        45       0.00029        95      1.32e-05      0.016      0.171       -0.00101    3.39e-06   83.2
        48       0.00029        99      1.32e-05      0.019      0.221       -0.00101    3.39e-06   83.2
        56       0.00029       207      1.32e-05      0.022      0.231       -0.00101    3.39e-06   83.2
        48       0.00029       140      5.22e-05      0.023      0.273       -0.00101    3.39e-06   83.2
        42       0.00029       121      5.22e-05      0.025      0.303       -0.00101    3.39e-06   83.2
        43       0.00029        78      1.32e-05      0.013      0.121       -0.00101    3.39e-06   83.2
        39         0.029       100      1.32e-05      0.015      0.151       -0.00101    3.39e-06   83.2
        41       0.00029       105      1.32e-05      0.016      0.171       -0.00101    3.39e-06   83.2
        43       0.00029       110      1.32e-05      0.019      0.221       -0.00101    3.39e-06   83.2
        45       0.00029       135      1.32e-05      0.022      0.231       -0.00101    3.39e-06   83.2
        55         0.228        70         6e-05      0.014       0.12         -0.002     1.8e-05   84.15
        60         0.228        75         5e-05      0.016       0.15         -0.002     1.8e-05   84.15
        50         0.228        85         6e-05      0.019       0.17         -0.002     1.8e-05   84.15
        50         0.228       110         7e-05       0.02        0.2         -0.002     1.8e-05   84.15
        58         0.228       100         7e-05      0.019       0.23         -0.003     1.8e-05   84.15
        66         0.228       105         7e-05      0.021       0.27         -0.005     1.5e-05   84.15
        70         0.228       100         7e-05      0.022        0.3         -0.005     1.5e-05   84.15
        34       0.00043       103       1.9e-05      0.015      0.151      -0.000614    7.78e-06   115.36
        31       0.00043       132      1.56e-05      0.017      0.201      -0.000295    7.78e-06   115.36
        300      0.000925      3.06       1.8e-05      0.046       0.36          0.021     9.8e-06  112.125
        100      0.000925      1.68       1.5e-05      0.087       0.25          0.004     9.8e-06  112.125
        85      0.000925      1.66         8e-06      0.058       0.36          0.008     9.8e-06   112.125
        263         0.002      33.5      2.12e-05      0.126       0.33        0.00419    9.76e-06  167.7
        165         0.002      18.8      2.96e-07       0.03       0.12         0.0002    6.65e-06  130.9111
        296         0.082      36.1      5.16e-05      0.171       0.33        0.00394    9.54e-06  130.9111];

  
    R1_L  = List_of_Inductors(Select_L,1);
    R2_L = List_of_Inductors(Select_L,2);
    C_L = List_of_Inductors(Select_L,3);
    k1 = List_of_Inductors(Select_L,4);
    k2  = List_of_Inductors(Select_L,5);
    k3 = List_of_Inductors(Select_L,6);
    k4 = List_of_Inductors(Select_L,7);
    k5 = List_of_Inductors(Select_L,8);
    L_area = List_of_Inductors(Select_L,9);

    C_L = C_L*1e-12;
    Rvar1 = k1*sqrt(fs);
    Rvar2 = k2*sqrt(fs);
    Lvar = (k3-k4*log(k5*fs))*1e-6;
    
    %penalty  = tot_area+L_area / 1;
    penalty = 0;
    %% This is all caluclations to set up the variables need to find the SS
    % Solution
 Vg = 20;
        Vbat1 = 4;
        Vbat2 = 4;
        M = (Vbat1+Vbat2)/(Vg);
        Vout = (Vbat1+Vbat2);
        Pout = 25;
        Vg = 20.5;
        Dboost = 0.15;
        Dbuck = M*(1-Dboost);
        Vout = M*Vg;
        Iout = Pout/Vout;
        
        Rg = 1e-3;
        Rbatt = 5e-3;
        Rcap = 2e-3;
        Cout = 2e-6;
        
        %{
        %L1 = 22e-6;
        
        %}
        
        %% Equation for Inductance
        % L1 has already been set by input
        %RL = L1*(9.32e-3/4.7e-6);
        
        
        Cfly = 23.5e-6;
        CflyESR = 2.2e-3;
        
        Vb1 = Vbat1;
        Vb2 = Vbat2;
        Rb = Rbatt;
        % EPC 2024
        %{
        r_on = 1.5e-3;
        Cds = 1620e-12;
        [M1_R,M2_R,M3_R,M4_R,M5_R,M6_R] = deal(r_on);
        [M1_C,M2_C,M3_C,M4_C,M5_C,M6_C] = deal(Cds);
        %}
        
        

        Ts = 1/fs;
        dt = 5e-9;
        
        u = [Vg Vbat1 Vbat2 1.5 1.5 1.5 1.5 1.5 1.5]';
            
        
        
    
    
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
    
    ts = [0.25 0.25 0.25 0.25]./Ts;
    
    % The inital guess of time intervals % The inital guess of time intervals
    % Assigned later dynamically
    
    
    % List all of the numerical components in the netlist file for all
    % FETs you must use the syntax used below:
    Numerical_Components = {
        'C1' Cfly
     'C2' Cout
     'C3' Cout
     'L1' Lvar
     'M1_C' Coss(1)
     'M2_C' Coss(2)
     'M3_C' Coss(3)
     'M4_C' Coss(4)
     'M5_C' Coss(5)
     'M6_C' Coss(6)
     'R1' R2_L
     'R2' Rbatt
     'R3' Rbatt
     'R4' Rcap
     'R5' Rcap
     'R6' CflyESR
     'R7' Rvar1
     'R8' R1_L
     'R9' Rvar2
     'C4' C_L
     'R10' Rg
     'M1_R' ron(1)
     'M2_R' ron(2)
     'M3_R' ron(3)
     'M4_R' ron(4)
     'M5_R' ron(5)
     'M6_R' ron(6)
     };
    
    % List out all char variables in the
    Switch_Resistors = {
        'M1_R'
        'M2_R'
        'M3_R'
        'M4_R'
        'M5_R'
        'M6_R'
        };
    % List of the switch sequency. Organized by: the FETs (column) vs time
    % interval (rows) matching Switch_Resistors and ts respectivly
    
    ON = 1;
    OFF = 0;
    Switch_Sequence = [
    ON  OFF ON  ON  OFF OFF   % Q1 ON
    OFF OFF ON  ON  OFF OFF  % Dead Q1 Q2
    OFF ON  ON  ON  OFF OFF  % Q2 ON
    OFF OFF ON  ON  OFF OFF  % Dead
    
    ];

    
    % List all the resistances of the diodes or FETS when they are on or
    % off
    SW_OFF = ones(1,6).*10000000;
    
    %SW_ON = [M1_R,M2_R,etc...]
    SW_ON = [ron(1) ron(2) ron(3) ron(4) ron(4) ron(3)];
    SW = [SW_OFF;SW_ON;SW_ON];
    
    Diode_Forward_Voltage = [1 1 1 1 1 1]'.*1.5;
    
    
    
    %% This creates and runs the parser to set up the converter
    
    parse = NetListParse();
    parse.initialize(filename,Voltage,Current);
    parse.Component_Values = Numerical_Components; %%% Added to test numerical stuff in code
    parse.cutset_loop_num();
    
    top = SMPStopology();
    top.Parser = parse;
    
    conv = SMPSconverter();
    conv.Topology = top;
    
    Order  = 1:4;
    
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
    sim.Converter = conv;
    sim_Rb = SMPSim;
    sim_Rb.Converter = conv;
    
    
    
    
    
    
    etaSim = [];
    PlossSim = [];
    PoutSim = [];
    conditions = [];
    PossSim = [];
    drange = linspace(0.39, .98, 12);
    Vscloc1 = 15;
    Vscloc2 = 18;
    Illoc = 1;
    Ib1loc = 2;
    Ib2loc = 3;
    VgPOS = 1;
    %h = waitbar(0, 'Simulating');
    Rg_place = 21;
    % Just for the buck boost not for the sc fib one 
    modSchemes = Switch_Sequence;
    
    for i = 1:size(modSchemes,3)
        %  for i = 3:6
        swvec = modSchemes(:,:,i);
        conv.Switch_Sequence = swvec;
        
        parse.Component_Values(Rg_place,2) = {1e6};
        conv.Element_Properties(Rg_place,2) = {1e6};
        
        sim_Rb=Run_SS_Converter_num_no_diode(sim_Rb,conv);
        
        parse.Component_Values(Rg_place,2) = {Rg};
        conv.Element_Properties(Rg_place,2) = {Rg};
        
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
        
        
        
        %  waitbar((i-1)/size(modSchemes,3), h);
        for d = drange
           % dt = 0.0025;
            ds = [d-dt1, dt1, 1-d-dt2, dt2];
            ts = ds/sum(ds)*Ts;
            conv.ts = ts;
            sim_Rb.ts = ts;
            sim.ts = ts;
            
            
            
           % Xss = sim_Rb.SS_Soln();
           % [ avgXs, avgYs ] = sim_Rb.ssAvgs(Xss);
            
            %OLVin = avgYs(Vscloc1)+avgYs(Vscloc2);
            OLVin = 8/d;
            % MaxVin = OLVin+1;
            MaxVin = OLVin+2;
            if OLVin<8 || OLVin>21
                continue
            end
            
            

            %{
        if i==2 && MaxVin>10
            continue
        end
        
        if i==3 && MaxVin>15
            continue
        end
        
        if i==4 && MaxVin>16 && OLVin<13
            continue
        end
        
        if i==5 && MaxVin>23 && OLVin<18
            continue
        end
        
        if i==6 && MaxVin>27 && OLVin<19
            continue
        end
            %}
            
            
            
            % waitbar((i-1)/size(modSchemes,3) + 1/size(modSchemes,3)*(find(drange==d)-1)/length(drange) , h);
            
            Vinrange = OLVin:.05:MaxVin;
            for Vin = Vinrange
                u(VgPOS) = Vin;
                conv.u = u;
                sim.u = u;
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
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
                
                if fail
                    
                    sim.As_saved = []  ;
                    sim.Bs_saved = [] ;
                    sim.Cs_saved = [] ;
                    sim.Ds_saved = [] ;
                    sim.u_saved = [];
                    sim.eigA_saved = [];
                    sim.ONorOFF_saved = [] ;
                    sim.ts_saved = [];
                    sim.Xs_saved = [];
                    
                end
                
                while fail && cycle < 2
                    
                    
                    
                    
                    fail = sim.Three_tier_diode_correct_num(iterations,0,1);
                    cycle = cycle+1;
                end
                
                %  if fail == true

               % end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                Xss = sim.SS_Soln();
                [ avgXs, avgYs ] = sim.ssAvgs(Xss);
                %{
            [ xs, t, ys ] = sim.SS_WF_Reconstruct;
            
            RMS_Iin = rms(ys(Illoc,:));
            RMS_Ib1 = rms(ys(Ib1loc,:));
            RMS_Ib2 = rms(ys(Ib2loc,:));
            
                
                
                
                %}
                
                %             avgYs(1:2)
                %     sim.plotAllStates(1)
                %     sim.plotAllOutputs(2)
                
                %             swActions = sum(abs(diff([swvec; swvec(1,:)],1)))/2;
                %             Vblock = [Vb1 Vb1 Vb1 Vb2 Vb2 Vb2 Vb1 Vb1 Vb2 Vb2 2*Vb1 2*Vb1 2*Vb1 2*Vb2 2*Vb2 2*Vb2];
                %             Poss = sum(1/2*Coss*Vblock.^2.*swActions*fs);
                
                
                %{
                 for the LM1205
                 peak source = 1.2 A
                peak sink = 5 A
                rise time = 7 ns
                fall time = 3.5 ns
                Qov = 1.4nC for 2055
                
                
 %}               
                
                Poverlap = abs(0.5*Vin*Xss(11,1)*1.4e-9/Ts)+abs(0.5*Vin*Xss(11,4)*1.4e-9/Ts);
                Poss = 0;
               % Poverlap = 0;
                % for when there is no rms calculation
                %%{
                eta = (avgYs(Ib1loc)*Vb1 + avgYs(Ib2loc)*Vb2 - Poss)/(Vin*-avgYs(Illoc) - avgYs(Ib1loc)^2*Rb - avgYs(Ib2loc)^2*Rb);
                Ploss = -(avgYs(Ib1loc)*Vb1 + avgYs(Ib2loc)*Vb2) + -(Vin*avgYs(Illoc)) + Poss+Poverlap;
                Pout = (avgYs(Ib1loc)*Vb1 + avgYs(Ib2loc)*Vb2) - Poss;
                %}
                
                % True input to output eff output of converter
                %{
            P_calc_Vin = Vin-(RMS_Iin*(RL*2));
            P_calc_Vout1 = Vb1+RMS_Ib1*Rb;
            P_calc_Vout2 = Vb2+RMS_Ib2*Rb;
            
            eta = (P_calc_Vout1*avgYs(Ib1loc) + P_calc_Vout2*avgYs(Ib2loc)) / -(P_calc_Vin*avgYs(Illoc));
            Ploss = -(P_calc_Vout1*avgYs(Ib1loc)+P_calc_Vout2*avgYs(Ib2loc))-(P_calc_Vin*avgYs(Illoc));
            Pout = (P_calc_Vout1*avgYs(Ib1loc)+P_calc_Vout2*avgYs(Ib2loc));
                %}
                if eta>1||eta<.3
                    continue
                end
                
                if Pout<0 || Ploss<0
                    continue
                end
                
                % for specific power
                %{
                if Pout<40 || Pout>60
                    continue
                end
                %}
                
                etaSim = [etaSim; eta];
                PlossSim = [PlossSim; Ploss];
                PoutSim = [PoutSim; Pout];
                %             PossSim = [PossSim; Poss];
                conditions = [conditions; d, i, Vin];
                
                if Pout > 70
                    break;
                end
                
                %             effSim(i,j) = eta;
                %             PlossSim(i,j) = Ploss;
                %             PoutSim(i,j) = (avgYs(Ib1loc)*Vb1 + avgYs(Ib2loc)*Vb2);
            end
            
            
            % %         %% Modeling
            % %         Ig = avgXs(Illoc);
            % %
            % %         dVc5 = Ig*sum(swvec(:,13-2).*ts')/Cfly2;
            % %         dVc1 = Ig*sum(swvec(:,1).*ts')/Cfly1;
            % %
            % %         Pdir = Vb1*Ig*sum(swvec(:,9-2).*ts')/Ts + Vb2*Ig*sum(swvec(:,12-2).*ts')/Ts;
            % %         Plossdir = 2*Ig^2*( ron*sum(swvec(:,9-2).*ts') + ...
            % %             ron*sum(swvec(:,13-2).*ts') + ...
            % %             2*ron*sum((1-swvec(:,13-2)).*ts') + ...
            % %             ron*sum(swvec(:,1).*ts') + ...
            % %             2*ron*sum((1-swvec(:,1)).*ts') + ...
            % %             Rb*Ts ...
            % %             )/Ts;
            % %
            % %
            % %         Psc = 2*1/2*Cfly2*((8+dVc5)^2 - 8^2)*fs + 2*1/2*Cfly1*((4+dVc1)^2 - 4^2)*fs;
            % %         PlossSc = 2*1/2*Cfly2*(dVc5^2)*fs + 2*1/2*Cfly1*(dVc1^2)*fs;
            % %
            % %         Poutsim = (avgYs(Ib1loc)*Vb1 + avgYs(Ib2loc)*Vb2);
            % %         Poutmodel = Pdir + Psc;
            % %
            % %     %     Plosssim = Ploss;
            % %     %     Plossmodel = Plossdir + PlossSc;
            % %
            % %         PoutModel(i,j) = Pdir + Psc;
            % %         PlossModel(i,j) =  Plossdir + PlossSc;
            % %         effModel(i,j) = Poutmodel/(Poutmodel + Plossdir + PlossSc);
            
            
        end
    end
    
    
    
    %close(h)
    
    % figure(4)
    % subplot(2,1,1)
    % plot(PoutSim, PlossSim, 'ok');
    % hold on;
    % plot(PoutModel, PlossModel)
    %
    % subplot(2,1,2)
    % plot(PoutSim, effSim, 'ok');
    % hold on;
    % plot(PoutModel, effModel)
    %
    
    % % [f,c] = contourf(repmat(Vgrange,length(drange),1),PoutSim,effSim*100);
    % % xlabel('Vg');
    % % ylabel('P_{out}');
    % % ylim([0 100])
    % % colorbar
    % % c.LevelList = [.8 .9 .95 .96 .97 .98 .99]*100;
    %{
mineff = .3;
etaSim(etaSim<mineff) = mineff;
etaSim(etaSim>1) = mineff;

locs = etaSim > mineff;

VgSim = conditions(locs,end);
F = scatteredInterpolant(VgSim, PoutSim(locs), etaSim(locs)*100, 'linear','none');

figure(101)

Vgrange = 5:.25:22;
PoutRange = 0:1:80;
[VgMesh, PoutMesh] = meshgrid(Vgrange, PoutRange);
[f,c] = contourf(Vgrange,PoutRange,F(VgMesh,PoutMesh),'ShowText','on');
xlabel('Vg');
ylabel('P_{out}');
% ylim([0 80])
colorbar
 c.LevelList = [.8 .9 .95 .96 .97 .98 .99]*100;

hold on;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mineff = .3;
etaSim(etaSim<mineff) = mineff;
etaSim(etaSim>1) = mineff;

locs = etaSim > mineff;

VgSim = conditions(locs,end);
F = scatteredInterpolant(VgSim, PoutSim(locs), PlossSim(locs), 'linear','none');

figure(104)

Vgrange = 5:.25:22;
PoutRange = 0:1:80;
[VgMesh, PoutMesh] = meshgrid(Vgrange, PoutRange);
[f,c] = contourf(Vgrange,PoutRange,F(VgMesh,PoutMesh),'ShowText','on');
xlabel('Vg');
ylabel('P_{out}');
% ylim([0 80])
colorbar
c.LevelList = [0 1 2 3 4 5 6 7 8 9];

hold on;


Vq = interp2(Vgrange,PoutRange,F(VgMesh,PoutMesh,,Yq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3)


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
colorbar
[c] = scatter(conditions(:,end),PoutSim,[],etaSim*100);
xlabel('Vg');
ylabel('P_{out}');
ylim([0 100])
colorbar

    %}
    
    
    
    
    
    mineff = .3;
    etaSim(etaSim<mineff) = mineff;
    etaSim(etaSim>1) = mineff;
    
    locs = etaSim > mineff;
    
    VgSim = conditions(locs,end);
    F = scatteredInterpolant(VgSim, PoutSim(locs), PlossSim(locs), 'linear','none');
    
    
    Vgrange = 10:.25:20;
    PoutRange = 20:1:60;
    
    sigma1 = 30;
    sigma2 = 750;
    
    mu = [16 50];
    Sigma = [sigma1, sqrt(sigma1*sigma2-1); sqrt(sigma1*sigma2-1), sigma2];
    
    [X1,X2] = meshgrid(Vgrange',PoutRange');
    X = [X1(:) X2(:)];
    
    p = mvncdf(X,mu,Sigma);
    
    Z = reshape(p,length(PoutRange),length(Vgrange));
    
    weighted_vals=F(X1,X2).*(Z+Z([length(PoutRange):-1:1],[length(Vgrange):-1:1]));
    
    
    Added_ploss = [];
    stick = [];
    stick = mean(mean(weighted_vals));
    stick = stick+penalty;
    
    
    if isnan(stick)||isinf(stick)||stick<0
        stick = 100;
    end
    
catch ME
    
    
    stick = 100;
    
end

return

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



