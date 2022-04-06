function [stick] = AURA_Eff_Sweep_IC_reduce_USB_SP_3(X)

%     _   _   _  ____    _
%    / \ | | | |/ _  |  / \
%   / _ \| | | | (_| | / _ \
%  / ___ | |_| |> _  |/ ___ \
% /_/   \_\___//_/ |_/_/   \_\



% Current konwn issues use FETs for diodes. (Just have them never turn
% on)

try
    
    % Place in the filename
    filename = 'SC_FIB_AURA_USB_SP.net'; % Place netlist filename here that you want to run
    
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
    
    
    
    
    penalty  = 0;
    adjust = [ones(1,length(X)-3).*20, 1e-6, 100, 100] ;
    X = X./adjust;
    
    W_selection = X(1:end-3);
    fs = X(end-2);
    dt1 = X(end-1);
    dt2 = X(end);
    ron = [];
    Coss = [];
    tot_area = 0;
    
    if length(W_selection)==8
        FET_selection = [2 2 2 3 3 3 3 3];
    elseif length(W_selection)==9
        FET_selection = [2 2 2 3 3 3 3 3 3];
    elseif length(W_selection)==5
        FET_selection = [3 3 3 3 3];
    else
        fprintf('Something went wrong /n')
    end
    
    
    
    
    %% Determining Coss and ron based on w
    
    for select_FET = 1:length(FET_selection)
        % Use typical value for Coss and max value for rds from data sheet
        % for inial look
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
            
        otherwise
            
            stick = 100;
            return
       
                
        end
        
    end
    
    
    
    %% This is all caluclations to set up the variables need to find the SS
    % Solution
    Vg = 12;
    
    
    Vb1 = 4;
    Vb2 = 4;
    
    VgPOS = 1;
    Vb1POS = 2;
    Vb2POS = 3;
    
    
    
    Rb = 5e-3;
    RL = 5e-3;
    
    
    Co = 20e-6; ESRo = 2e-3;
    %Cfly1 = 5e-6; ESR1 = 3e-3;  % 5V Cap
    %Cfly2 = 2.5e-6; ESR2 = 3e-3; % 10V Cap
    Cfly1 = 10e-6; ESR1 = 3e-3;  % 5V Cap
    Cfly2 = 10e-6; ESR2 = 3e-3; % 10V Cap
    
    Lc = 100e-9;
    Ts = 1/fs;
    
    
    RL1 = 35;
    LL1 = 300e-9;
    CL1 = 375e-11;
    
    
    RL2 = 0.01;
    LL2 = 900e-9;
    
    
    
    
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
        0     1     1     1     0     1     1     0     0     1     0     1     1     1     0     1
        0     1     1     1     0     1     1     0     0     1     0     0     0     1     0     1
        0     1     1     1     0     1     1     0     0     1     1     0     0     1     0     1
        0     1     1     1     0     1     1     0     0     1     0     0     0     1     0     1
        0     1     1     1     0     1     1     0     0     1     0     1     1     1     0     1
        0     1     1     1     0     1     1     0     0     1     0     1     1     0     0     0
        0     1     1     1     0     1     1     0     0     1     0     1     1     0     1     0
        0     1     1     1     0     1     1     0     0     1     0     1     1     0     0     0
        ];
    
    modSchemes(:,:,4) = [
        0     1     1     1     0     1     1     0     0     1     1     0     0     1     0     1
        0     0     0     1     0     1     1     0     0     1     1     0     0     1     0     1
        1     0     0     1     0     1     1     0     0     1     1     0     0     1     0     1
        0     0     0     1     0     1     1     0     0     1     0     0     0     0     0     0
        0     1     1     1     0     1     1     0     0     1     0     1     1     0     1     0
        0     1     1     0     0     0     1     0     0     1     0     1     1     0     1     0
        0     1     1     0     1     0     1     0     0     1     0     1     1     0     1     0
        0     1     1     0     0     0     1     0     0     1     0     0     0     0     0     0
        ];
    
    modSchemes(:,:,5) = [
        1     0     0     1     0     1     1     0     0     1     1     0     0     1     0     1
        1     0     0     1     0     1     1     0     0     1     1     0     0     0     0     0
        1     0     0     1     0     1     1     0     0     1     1     0     0     0     1     0
        0     0     0     0     0     0     1     0     0     1     0     0     0     0     1     0
        0     1     1     0     1     0     1     0     0     1     0     1     1     0     1     0
        0     1     1     0     1     0     1     0     0     1     0     0     0     0     1     0
        0     1     1     0     1     0     1     0     0     1     1     0     0     0     1     0
        0     0     0     0     0     0     1     0     0     1     1     0     0     0     0     0
        ]; 
    
     modSchemes(:,:,6) = [
        1     0     0     1     0     1     1     0     0     1     1     0     0     1     0     1
        1     0     0     0     0     0     1     0     0     1     1     0     0     0     0     0
        1     0     0     0     1     0     1     0     0     1     1     0     0     0     1     0
        0     0     0     0     1     0     1     0     0     1     0     0     0     0     1     0
        0     1     1     0     1     0     1     0     0     1     0     1     1     0     1     0
        0     0     0     0     1     0     1     0     0     1     0     0     0     0     1     0
        1     0     0     0     1     0     1     0     0     1     1     0     0     0     1     0
        1     0     0     0     0     0     1     0     0     1     1     0     0     0     0     0
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
    
    ts = [0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25]./Ts;
    
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
        'L1' LL2
        'L2' LL1
        'R3' RL2
        'R4' RL1
        'C7' CL1
        'R1' Rb
        'R2' Rb
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
        'M11_C' Coss(1)
        'M12_C' Coss(2)
        'M13_C' Coss(3)
        'M14_C' Coss(2)
        'M15_C' Coss(1)
        'M16_C' Coss(3)
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
        'M11_R' ron(1)
        'M12_R' ron(2)
        'M13_R' ron(3)
        'M14_R' ron(2)
        'M15_R' ron(1)
        'M16_R' ron(3)
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
    SW_ON = [ron(1) ron(2) ron(3) ron(2) ron(1) ron(3) ron(4) ron(5) ron(5) ron(4) ron(1) ron(2) ron(3) ron(2) ron(1) ron(3)];
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
    
    Order  = 1:8;
    
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
    drange = linspace(0.05, .95, 8);
    Vscloc1 = 39;
    Vscloc2 = 40;
    Illoc = 1;
    Ib1loc = 2;
    Ib2loc = 3;
    %h = waitbar(0, 'Simulating');
    
    for i = 1:size(modSchemes,3)
        %  for i = 3:6
        swvec = modSchemes(:,:,i);
        conv.Switch_Sequence = swvec;
        
        parse.Component_Values(9,2) = {1e6};
        conv.Element_Properties(9,2) = {1e6};
        
        sim_Rb=Run_SS_Converter_num_no_diode(sim_Rb,conv);
        
        parse.Component_Values(9,2) = {RL2};
        conv.Element_Properties(9,2) = {RL2};
        
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
        
        if i==4
            J = 45645;
        end
        
        %  waitbar((i-1)/size(modSchemes,3), h);
        for d = drange
            
            %if d>0.5 && i==6
            %    continue
            %end
            if d<0.1 && i==5
                continue
            end
            
            if i==5
                continue
            end
                
           % dt = 0.0025;
            ds = [d-dt1, dt1, 1-d-dt2, dt2, d-dt1, dt1, 1-d-dt2, dt2];
            ts = ds/sum(ds)*Ts;
            conv.ts = ts;
            sim_Rb.ts = ts;
            sim.ts = ts;
            
            
            
            Xss = sim_Rb.SS_Soln();
            [ avgXs, avgYs ] = sim_Rb.ssAvgs(Xss);
            
            OLVin = avgYs(Vscloc1)+avgYs(Vscloc2);
            % MaxVin = OLVin+1;
            MaxVin = OLVin+2;
            if OLVin<4 || OLVin>23
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
            
            Vinrange = OLVin:.1:MaxVin;
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
                iterations = 20;
                cycle = 0;
                fail = 1;
                %fprintf('--------------\n')
              %  fail = sim.Three_tier_diode_correct_num(iterations,0,0);
              %  
              %  while fail && cycle < 1
              %      fail = sim.Three_tier_diode_correct_num(iterations,0,1);
              %      cycle = cycle+1;
              %      
              %  end
                
               
                
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
                Poss = 0;
                
                % for when there is no rms calculation
                %%{
                eta = (avgYs(Ib1loc)*Vb1 + avgYs(Ib2loc)*Vb2 - Poss)/(Vin*-avgYs(Illoc) - avgYs(Ib1loc)^2*Rb - avgYs(Ib2loc)^2*Rb);
                Ploss = -(avgYs(Ib1loc)*Vb1 + avgYs(Ib2loc)*Vb2) + -(Vin*avgYs(Illoc)) + Poss;
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
                
                if Pout > 50
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

figure(103)

Vgrange = 5:.25:22;
PoutRange = 0:1:80;
[VgMesh, PoutMesh] = meshgrid(Vgrange, PoutRange);
[f,c] = contourf(Vgrange,PoutRange,F(VgMesh,PoutMesh),'ShowText','on');
xlabel('Vg');
ylabel('P_{out}');
% ylim([0 80])
colorbar
c.LevelList = [0 0.25 0.5 1 1.5 2 3 4 5];

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
    
    
    Vgrange = 7:.25:18;
    PoutRange = 10:1:40;
    
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
    stick = mean(weighted_vals,'all');
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



