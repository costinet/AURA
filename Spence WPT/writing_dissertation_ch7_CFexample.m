%% ------------------------------------------------------------- 
% -------------- 7 Level Model Anchored at Vrec ---------------- 
% -------------------------------------------------------------- 
clear all; 
close all; 
clc; 
tic;
addpath('functions');
addpath('dataSave');
addpath('includes');
addpath('simModels');
fntSz = 14;
opts = bodeoptions('cstprefs');     opts.FreqUnits = 'kHz'; 
    opts.XLabel.FontSize = fntSz;   opts.YLabel.FontSize = fntSz; 
    opts.Title.FontSize = fntSz;    opts.TickLabel.FontSize = 12; 
    opts.PhaseWrapping = 'off';                 % wrap from -180 to 180
    opts.PhaseMatching = 'on';                  % set phase +- 360 deg from ... 
    opts.PhaseMatchingFreq = 0;                 % "set" phase at this freq
    opts.PhaseMatchingValue = -91;              % set phase at [^freq] near to this val


% ------------- Save System Characteristic Data ----------------
if(1)
% Tank Tuning, June ~15th, 2019
Ltx = 10.78e-6;                     Lrx = 12.11e-6; 
k = 0.5;                            Lm = k*sqrt(Ltx*Lrx);
Ls = Lrx - Lm;                      Lp = Ltx - Lm; 
Cp = 221.0e-9;                      Cs = 95.13e-9;
Rp = 79.69e-3 + 75e-3;              Rs = 390.16e-3;
% Rp = 250e-3;                        Rs = 500e-3; 

% % CHONGWENS TANK
% Ltx = 10e-6; 
% Lrx = 10e-6; 
% k = 0.7; 
% Lm = k*sqrt(Ltx*Lrx); 
% Lp = Ltx - Lm; 
% Ls = Lrx - Lm;  
% Cp = 375e-9; 
% Cs = 375e-9; 
% Rp = 200e-3; 
% Rs = 200e-3; 

% % ORIGINAL TANK (for testing)
% % Vg = 15.9931; %15.66
% % Vs = 14.6714;%5.01*3.03;
% Cp = 186.43e-9;%220.8e-9; %100.8e-9; 
% Cs = 174.09e-9;%179.2e-9;  %100.2e-9;  %179.2e-9;
% Rp = .1206;%66.9e-3 + 75e-3; 
% Rs = .61649;%348.9e-3; %70.9e-3;
% Lm = 7.4735e-6;%k*sqrt(Ltx*Lrx); 
% Lp = 6.4153e-6;%Ltx - Lm;
% Ls = 6.17365e-6;%Lrx - Lm;
% Rload = 18.5;                     % EDIT VALUES ON LOWER LINES!! 
% Vinv = 16.5; 
end
% --------------------------------------------------------------


% ---------------------- Setup Simulaiton ----------------------
if(1)
    w = 2*pi*150e3; 
    
    % Wednesday Jan 22nd (measured with coils taped together) 
    Ltx = 13.22e-6;                         Lrx = 13.26e-6;
    k = 0.773;                              Lm = k*sqrt(Ltx*Lrx);
    Ls = Lrx - Lm;                          Lp = Ltx - Lm;
    Cp = 757.56e-9;                         Cs = 111.76e-9;
%     Cp = 1/w^2/Ltx;                         Cs = 1/w^2/Lrx;
    Rp = 60.94e-3 + 8.20e-3 + 2*7.7e-3;  	Rs = (63.41e-3 + 41.39e-3);     % Bulk parasitic.
%     Rp = 60.94e-3 + 8.20e-3;  	Rs = 63.41e-3;                  % PCB parasitics. 

    if(1) % '1' if Rload, '0' if Iload in PLECS.
        % [Max Rload,    Min Rload,       Rload_0, itterate R, ---]
        RlMax = 100; RlMin = 0.025; Rload = RlMin; fRoutg = 1; fIoutg = 0; 
%         RlMax = 100; RlMin = 0.025; Rload = 001.5; fRoutg = 0; fIoutg = 0;     
    else
        % [Max Iload,   Min Iload,         Iload_0, itterate I, ---]
        IlMax = 005; IlMin = 0.01; Iload = IlMax/2; fIoutg = 1; fRoutg = 0; end	
    
%         Rload = 1.8; fRoutg = 0; fIoutg = 0; 

    
    indVout = 7;    indIp = 10;     indIs = 11;     indCp = 8;      indCs = 9; 
%     indVout = 7;     indIp = 3;      indIs = 4;     indCp = 1;      indCs = 2; 
    
    % Inverter         goal Vout         Vout tolerance
    Vinv = 10;      Voutg = 2.0;        Voutt = 0.0001;    
    Vinv = 15;      Voutg = 5.0;        Voutt = 0.005;           
    U = Vinv;                       	% Rload.
%     U = [Vinv, Iload]';           	% Iload.
%     U = [Vinv, repmat(Voutg,1,9)]';     % V src Flying Caps.

    sc2 = 0;        	% Figure print offset ('1920' for screen #2).

    ph = 150;           % Phase Vs-Vp.
    b11Mod = 811;       % Modulation Index.
    ph = 130;           % Phase Vs-Vp.
    b11Mod = 620;       % Modulation Index.

    ns = 11;         	% [11, -,11,37,13,28] 	# of states in system. 
    meters = 3;     	% [ 3, -, 2, 6, 3, 6] 	# number of outputs (PLECS meters). 
    kint = 16;         	% [16,16,16,16,16,16]  	# of intervals in system. 

    %         [  Rload,     ---,   Iload,  Par:R+L,   Par:R,  Vsrc FlyCaps]
    circNum = [   {''},  {'_2'},  {'_3'},   {'_4'},  {'_5'},        {'_6'}];
    circNum = circNum{1};      
    myPath = ['PLECS_7lvl/Circuit' circNum];      exString = ['Circuit' circNum]; 
end
% --------------------------------------------------------------         


% ------------------ Simulation Parameters ---------------------  
if(1)
    f = 150e3;            	
    Ts = 1/f;            	
    clks = 1000;                    % clocks / fundamental period.
    clk2 = clks/2; 
    Pclk150 = Ts/clks;              % Period of 150 MHz clock.
    fclk150 = 1/Pclk150;            % Frequency of 150 MHz clock.
    fclk20 = 20e6; 
    Pclk20 = 1/fclk20; 
    T_simulation = ...
        1/lcm(fclk20,fclk150);   

    nm = 3;                         % [#]       Levels in phase single leg.
    Vgs = 5;                        % [V]       Driving voltage. 
    Coss = 525e-12;                 % [F]       (datasheet) (Rough avg. from 0V-to-5V Vds)
        Qoss = Coss*Voutg;          % [C]       Qoss given delatV & ~Coss
    Qgs = 1.4e-9;                   % [C]       (datasheet) guess @ 10V (Qgs=1.5nC @ Vgs=12V)
    Igd = 1;                        % [A]       1 Ohm gate res. @ Vgs = 10V. 
    ton = Qgs/Igd;                  % [sec]     Overlap loss time.
        toff = ton; 
    Qrr = 5e-9;                     % [C]       (datasheet) *defined by design, not test*, If = 9A
        trr = ton;                          
    Vd = 0.86; 
    Qgate = 8.0e-9;                 % [C]       (datasheet) Qgate @ Vgs = 10
    Ron = 7.7e-3;                   % [Ohm]     (datasheet) FET on Res, @ Vgs = 10
%     Ron = 10e-3;                   % [Ohm]     (datasheet) FET on Res, @ Vgs = 10
    Resr = 1.39e-3;               	% [Ohm]     ESR Power Caps, Impedance Analyzer: 2020-01-24
%     Resr = 5e-3;               	% [Ohm]     ESR Power Caps, Impedance Analyzer: 2020-01-24
    Resr_out = Resr;                % [Ohm]     ESR at of output capacitance
    Rpar = 200e1;                  	% [Ohm]     PLECs parasitic resistance. 
    Cpower = 15.66e-6;            	% [F]     	Impedance Analyzer: 2020-01-24
    Cpower = 10e-6
      CA1 = Cpower; CA2 = Cpower; 
      CB1 = Cpower; CB2 = Cpower;
      Cout_near = Cpower; 
      Cout_far = Cpower *1.75; 
      Cout = Cout_far + Cout_near*2;    % (3.75 * Cpower)
      Ceq3 = Cpower*2 + Cout; 
      Cout = 50e-6; 
      
    dRon = 0e-3; 
    Ron_1H = Ron + 2*dRon;
    Ron_1L = Ron + 1*dRon;
    Ron_2H = Ron + 2*dRon;
    Ron_2L = Ron + 1*dRon;
    Ron_3H = Ron + 2*dRon;
    Ron_3L = Ron + 1*dRon;
    Ron_4H = Ron + 4*dRon;
    Ron_4L = Ron + 3*dRon;

    fact = 1;                                           fact2 = 1; 
    Rp1  = 0.43e-3 *fact;   	Lp1  = 1.00e-9;         RpL1  = 10000 *fact2; 
    Rp3  = 1.68e-3 *fact;       Lp3  = 3.16e-9;         RpL3  = 10000 *fact2;
    Rp4  = 0.45e-3 *fact;    	Lp4  = 0.98e-9;         RpL4  = 10000 *fact2; 
    Rp6  = 0.94e-3 *fact;     	Lp6  = 2.19e-9;         RpL6  = 10000 *fact2; 
    Rp7  = 0.44e-3 *fact;      	Lp7  = 0.99e-9;         RpL7  = 10000 *fact2; 
    Rp9  = 2.50e-3 *fact;      	Lp9  = 4.47e-9;         RpL9  = 10000 *fact2; 
    Rp10 = 0.52e-3 *fact;    	Lp10 = 1.17e-9;         RpL10 = 10000 *fact2; 
    Rp12 = 0.78e-3 *fact;    	Lp12 = 1.89e-9;         RpL12 = 10000 *fact2; 
    Rp13 = 0.38e-3 *fact;    	Lp13 = 0.90e-9;         RpL13 = 10000 *fact2; 
    Rp15 = 0.87e-3 *fact;     	Lp15 = 2.14e-9;         RpL15 = 10000 *fact2; 
    Rp16 = 0.55e-3 *fact;     	Lp16 = 1.24e-9;         RpL16 = 10000 *fact2; 
    Rp18 = 0.14e-3 *fact;     	Lp18 = 0.32e-9;         RpL18 = 10000 *fact2; 
    Rp19 = 3.01e-3 *fact;     	Lp19 = 14.94e-9;    	RpL19 = 10000 *fact2; 
    Rp20 = 2.57e-3 *fact;     	Lp20 = 8.85e-9;         RpL20 = 10000 *fact2; 
    Rp22 = 4.95e-3 *fact;      	Lp22 = 16.11e-9;      	RpL22 = 10000 *fact2; 
    Rp25 = 1.10e-3 *fact;    	Lp25 = 14.26e-9;        RpL25 = 10000 *fact2; 
    Rpj  = 0.26e-3 *fact;      	Lpj  = 5.03e-9;         RpLj  = 10000 *fact2; 

    Rtest = 100e-3; 

    % b12Mod = 620;
    % ph = 73; 

    % ph = 120;                       % phase shift (Control: Vp-Vs) [clks]
    % b11Mod = 900;

    % ph = 100;                       % phase shift (Control: Vp-Vs) [clks]
    % b11Mod = 620; 

    % ph = 150; 
    % b11Mod = 811; 
    % b11Mod = 820; 

    [M1,M2,M3,nVs,duty] = b12toMod(b11Mod);
    % [M1,M2,M3] = fundToMs(4/pi*5);
    % M1 = 19; 
    % M2 = 59; 
    % M3 = 107;
end
% --------------------------------------------------------------


% -------------- Extract State Spaces from PLECS ---------------
if(1)
    ssOrder = plecs('get', myPath, 'StateSpaceOrder');

    % [[columns = switch, rows = interval]]
    % "switches" below should match command: "ssOrder.Switches"
    % switches: [SA1H SA1L SA2H SA2L SA3H SA3L SA4H SA4L ...
    %           ... SB1H SB1L SB2sH SB2L SB3H SB3L SB4H SB4L]
    % intervals: [0 1A 2A 3A 1B 2B 3B]'
    order = [ 0 1 0 1 0 1 1 1   0 1 0 1 0 1 1 1;    % [0]
              0 1 0 1 1 0 0 0   0 1 0 1 0 1 1 1;    % [1A]
              0 1 1 0 1 0 0 0   0 1 0 1 0 1 1 1;    % [2A]
              1 0 1 0 1 0 0 0   0 1 0 1 0 1 1 1;    % [3A]
              0 1 0 1 0 1 1 1   0 1 0 1 1 0 0 0;    % [1B]
              0 1 0 1 0 1 1 1   0 1 1 0 1 0 0 0;    % [2B]
              0 1 0 1 0 1 1 1   1 0 1 0 1 0 0 0; ];	% [3B]
%     order = [ 0 1 0 1 0 1 0 0   0 1 0 1 0 1 0 0;    % [0]
%               0 1 0 1 1 0 0 0   0 1 0 1 0 1 1 1;    % [1A]
%               0 1 1 0 1 0 0 0   0 1 0 1 0 1 1 1;    % [2A]
%               1 0 1 0 1 0 0 0   0 1 0 1 0 1 1 1;    % [3A]
%               0 1 0 1 0 1 1 1   0 1 0 1 1 0 0 0;    % [1B]
%               0 1 0 1 0 1 1 1   0 1 1 0 1 0 0 0;    % [2B]
%               0 1 0 1 0 1 1 1   1 0 1 0 1 0 0 0; ];	% [3B]
%     order = [ 0 1 0 1 0 1 0 1   0 1 0 1 0 1 0 1;    % [0]           % High side turns on LATER. 
%               0 1 0 1 1 0 0 0   0 1 0 1 0 1 1 1;    % [1A]
%               0 1 1 0 1 0 0 0   0 1 0 1 0 1 1 1;    % [2A]
%               1 0 1 0 1 0 0 0   0 1 0 1 0 1 1 1;    % [3A]
%               0 1 0 1 0 1 1 1   0 1 0 1 1 0 0 0;    % [1B]
%               0 1 0 1 0 1 1 1   0 1 1 0 1 0 0 0;    % [2B]
%               0 1 0 1 0 1 1 1   1 0 1 0 1 0 0 0; ];	% [3B]
    top = topology;
    top.loadPLECsModel(myPath,order);
    save('PLECs_stateSpace.mat'); 
    % load('PLECs_stateSpace.mat'); 
end
% --------------------------------------------------------------


% ---------- Printing State Space Order -----------
if(1)
    disp('------------------------------------'); 
    fprintf('x =                U =            y = \n'); 
    for i = 1:1:length(ssOrder.States)
        if(contains(ssOrder.States{i},[exString '/WPT Coils:']))
            fprintf('  [%2.0f]  L%-10s', i, erase(ssOrder.States{i}, [exString '/WPT Coils:'])); 
        else
            fprintf('  [%2.0f]  %-11s', i, erase(ssOrder.States{i}, [exString '/'])); end
        
        if(i<=length(ssOrder.Inputs))
            fprintf('  [%.0f]  %-7s', i, erase(ssOrder.Inputs{i}, [exString '/'])); 
        else
            fprintf('              '); end
        
        if(i<=length(ssOrder.Outputs))
            fprintf('   [%.0f]  %-7s', i, erase(ssOrder.Outputs{i}, [exString '/'])); end
        fprintf('\n'); 
    end
    disp('------------------------------------'); 
    disp(['Indices used for calculations:']); 
    disp('Vout  Vcp  Vcs   ip   is'); 
    fprintf('%4.0f %4.0f %4.0f %4.0f %4.0f\n', indVout, indCp, indCs, indIp, indIs); 
    disp('------------------------------------'); 
end
% -------------------------------------------------


% ------------- Find Iload Such That Iout = Ioutg --------------
if(fIoutg) 
    % --------------- Search: Change Rload, Find Voutg -------------
    % Assumes monotonic. Has safeties for run-away. 
    itts = 0; ittsMax = 100; StpSz = 0.5; dirTog = -1;   	% Initialize itteration.
    err = 1e3; serr = 1; search = 1;  
    while(search)
        
        % ------------------ Save, Check, Increment --------------------
        sePrev = serr; ePrev = err;     % Save previous iteration.
        IlPrev = Iload; 
        Iload = Iload + StpSz*dirTog;  	% Step Iload appropriately.
        if(Iload>IlMax) Iload = IlMax; 
            elseif(Iload<IlMin) Iload = IlMin; end
        itts = itts + 1; 
%         U = [Vinv; Iload]; 
        U(2) = Iload; 
        % --------------------------------------------------------------

        % ------------------ Evaluate the Topology ---------------------
        top.loadPLECsModel(myPath,order);                     	% Solve new state space.    
        [td,te,te2,topA,topB,inPair,inPolD,inPolP,Vrs,Us] = ... % Solve modulation indices 
            modInts7(clk2,Pclk150,ph,M1,M2,M3,top.As(:,:,1),... % ... for this operating 
            top.As(:,:,2), top.As(:,:,3), top.As(:,:,4), ...  	% ... point (mod + ph). 
            top.As(:,:,5), top.As(:,:,6), top.As(:,:,7), ...
            top.Bs(:,:,1), top.Bs(:,:,2), top.Bs(:,:,3), ...
            top.Bs(:,:,4), top.Bs(:,:,5), top.Bs(:,:,6), ...
            top.Bs(:,:,7)); 
        Xss_all = stSt_solve(topA, topB, td, U);                    % Steady state solution. ('U' must be constant)
        Xss = Xss_all(:,1);                                         % Isolate beginning of period.
        % --------------------------------------------------------------

        % ------------------- Calc. Operating Point --------------------
        Xact = zeros(size(squeeze(topA(1,:,:)))); 
        for k = find(td/Pclk150~=0)                         	% Each interval.
            Xact(:,k) = 1/td(k) * (topA(:,:,k)\((expm(...       % Average of 'x' w/in
              topA(:,:,k)*td(k))-eye(ns))*Xss_all(:,k) +...     % ... interval 'k'. 
              topA(:,:,k)\(expm(topA(:,:,k)*td(k))-eye(ns)...
              -topA(:,:,k)*td(k))*topB(:,:,k)*U));
        end 
        Vout = sum(Xact(indVout,:).*td/Ts);    
        if(isnan(Vout)) error('Vout is NaN!!!'); end
        err = Voutg - Vout;                                     % Error.
        serr = sign(err);                                       % Sign of error.   
        % --------------------------------------------------------------
        
        % --------------------- Search Decisions -----------------------
        if(abs(err)<Voutt)                          % SOL'N FOUND.
            search = 0; valid = 1; msg = "[]";
        elseif(serr~=sePrev)                        % REVERSE DIRECTION.
            StpSz = StpSz/2; dirTog = -1*dirTog; 
        elseif(Iload==IlMax)                        % END REACHED. (no sol'n)
            search = 0; valid = 0; msg = "Max";
        elseif(Iload==IlMin && itts~=1)             % BEGINNING REACHED. (no sol'n)
            search = 0; valid = 0; msg = "Min";
        elseif(itts > ittsMax)                      % TOO MANY ITTERATIONS.
            search = 0; valid = 0; msg = "its"; end
        % --------------------------------------------------------------            
    end
    % -------------------------------------------------------------- 
end
% --------------------------------------------------------------


% ------------- Find Rload Such That Vout = Voutg --------------
if(fRoutg) 
    % --------------- Search: Change Rload, Find Voutg -------------
    % Assumes monotonic. Has safeties for run-away. 
    itts = 0; ittsMax = 100; StpSz = 10; dirTog = 1;   	% Initialize itteration.
    err = 1e3; serr = 1; search = 1; Rload = RlMin; 
    while(search)

        % ------------------ Save, Check, Increment --------------------
        sePrev = serr; ePrev = err;     % Save previous iteration.
        RlPrev = Rload; 
        Rload = Rload + StpSz*dirTog;  	% Step Rload appropriately.
        if(Rload>RlMax) Rload = RlMax; 
            elseif(Rload<RlMin) Rload = RlMin; end
        itts = itts + 1; 
        % --------------------------------------------------------------

        % ------------------ Evaluate the Topology ---------------------
        top.loadPLECsModel(myPath,order);                     	% Solve new state space.    
%         [td,te,te2,topA,topB,inPair,inPolD,inPolP,Vrs,Us] = ... % Solve modulation indices 
%             modInts7(clk2,Pclk150,ph,M1,M2,M3,top.As(:,:,1),... % ... for this operating 
%             top.As(:,:,2),top.As(:,:,3),top.As(:,:,4),...       % ... point (mod + ph). 
%             top.As(:,:,5),top.As(:,:,6),...
%             top.As(:,:,7),top.Bs(:,:,1)); 
     	[td,te,te2,topA,topB,inPair,inPolD,inPolP,Vrs,Us] = ... % Solve modulation indices 
        modInts7(clk2,Pclk150,ph,M1,M2,M3,top.As(:,:,1),... % ... for this operating 
        top.As(:,:,2), top.As(:,:,3), top.As(:,:,4), ...  	% ... point (mod + ph). 
        top.As(:,:,5), top.As(:,:,6), top.As(:,:,7), ...
        top.Bs(:,:,1), top.Bs(:,:,2), top.Bs(:,:,3), ...
        top.Bs(:,:,4), top.Bs(:,:,5), top.Bs(:,:,6), ...
        top.Bs(:,:,7));   
        Xss_all = stSt_solve(topA, topB, td, U);                    % Steady state solution. ('U' must be constant)
        Xss = Xss_all(:,1);                                         % Isolate beginning of period.
        % --------------------------------------------------------------

        % ------------------- Calc. Operating Point --------------------
        Xact = zeros(size(squeeze(topA(1,:,:)))); 
        for k = find(td/Pclk150~=0)                           	% Each interval.
            Xact(:,k) = 1/td(k) * (topA(:,:,k)\((expm(...       % Average of 'x' w/in
              topA(:,:,k)*td(k))-eye(ns))*Xss_all(:,k) +...     % ... interval 'k'. 
              topA(:,:,k)\(expm(topA(:,:,k)*td(k))-eye(ns)...
              -topA(:,:,k)*td(k))*topB(:,:,k)*U)); end 
        Vout = sum(Xact(indVout,:).*td/Ts);    
        if(isnan(Vout)) error('Vout is NaN!!!'); end
        err = Voutg - Vout;                     % Error.
        serr = sign(err);                       % Sign of error.   
        % --------------------------------------------------------------
        
        % --------------------- Search Decisions -----------------------
        if(abs(err)<Voutt)                          % SOL'N FOUND.
            search = 0; valid = 1; msg = "[]";
        elseif(serr~=sePrev)                        % REVERSE DIRECTION.
            StpSz = StpSz/2; dirTog = -1*dirTog; 
        elseif(Rload==RlMax)                        % END REACHED. (no sol'n)
            search = 0; valid = 0; msg = "Max";
        elseif(Rload==RlMin && itts~=1)             % BEGINNING REACHED. (no sol'n)
            search = 0; valid = 0; msg = "Min";
        elseif(itts > ittsMax)                      % TOO MANY ITTERATIONS.
            search = 0; valid = 0; msg = "its"; end
        % --------------------------------------------------------------            
    end
    % -------------------------------------------------------------- 
end
% --------------------------------------------------------------


% ------------------ Evaluate the Topology ---------------------
if(1)
    top.loadPLECsModel(myPath,order);                     	% Solve new state space. 
%     [td,te,te2,topA,topB,inPair,inPolD,inPolP,Vrs,Us] = ... % Solve modulation indices 
%         modInts7(clk2,Pclk150,ph,M1,M2,M3,top.As(:,:,1),... % ... for this operating 
%         top.As(:,:,2),top.As(:,:,3),top.As(:,:,4),...       % ... point (mod + ph). 
%         top.As(:,:,5),top.As(:,:,6),...
%         top.As(:,:,7),top.Bs(:,:,1)); 
    [td,te,te2,topA,topB,inPair,inPolD,inPolP,Vrs,Us] = ... % Solve modulation indices 
        modInts7(clk2,Pclk150,ph,M1,M2,M3,top.As(:,:,1),... % ... for this operating 
        top.As(:,:,2), top.As(:,:,3), top.As(:,:,4), ...  	% ... point (mod + ph). 
        top.As(:,:,5), top.As(:,:,6), top.As(:,:,7), ...
        top.Bs(:,:,1), top.Bs(:,:,2), top.Bs(:,:,3), ...
        top.Bs(:,:,4), top.Bs(:,:,5), top.Bs(:,:,6), ...
        top.Bs(:,:,7));   
    Xss_all = stSt_solve(topA, topB, td, U);                % Steady state solution. ('U' must be constant)
    Xss = Xss_all(:,1);                                     % Isolate beginning of period.
    for k = find(td/Pclk150~=0)
        Xact(:,k) = 1/td(k) * (topA(:,:,k)\((expm(...     	% Average of 'x' within 
          topA(:,:,k)*td(k))-eye(ns))*Xss_all(:,k) +...    	% ... interval 'k'.
          topA(:,:,k)\(expm(topA(:,:,k)*td(k))-...        	% ------------------------------
          eye(ns)-topA(:,:,k)*td(k))*topB(:,:,k)*U)); end  
    Vout = sum(Xact(indVout,:).*td/Ts);                     % Period average of Vout.
    
    Vcb3 = sum(Xact(1,:).*td/Ts); 
    Vca1 = sum(Xact(2,:).*td/Ts); 
    Vca2 = sum(Xact(3,:).*td/Ts); 
    Vca3 = sum(Xact(4,:).*td/Ts); 
    Vcb1 = sum(Xact(5,:).*td/Ts); 
    Vcb2 = sum(Xact(6,:).*td/Ts); 
            
    topC = zeros(meters,ns,k); 
        inds = []; inds = find(Vrs==0);  for i = 1:1:length(inds) topC(:,:,inds(i)) = top.Cs(:,:,1); end
        inds = []; inds = find(Vrs==1);  for i = 1:1:length(inds) topC(:,:,inds(i)) = top.Cs(:,:,2); end
        inds = []; inds = find(Vrs==2);  for i = 1:1:length(inds) topC(:,:,inds(i)) = top.Cs(:,:,3); end
        inds = []; inds = find(Vrs==3);  for i = 1:1:length(inds) topC(:,:,inds(i)) = top.Cs(:,:,4); end
        inds = []; inds = find(Vrs==-1); for i = 1:1:length(inds) topC(:,:,inds(i)) = top.Cs(:,:,5); end
        inds = []; inds = find(Vrs==-2); for i = 1:1:length(inds) topC(:,:,inds(i)) = top.Cs(:,:,6); end
        inds = []; inds = find(Vrs==-3); for i = 1:1:length(inds) topC(:,:,inds(i)) = top.Cs(:,:,7); end
        
    topD = zeros(meters,length(U),k); 
        inds = []; inds = find(Vrs==0);  for i = 1:1:length(inds) topD(:,:,inds(i)) = top.Ds(:,:,1); end
        inds = []; inds = find(Vrs==1);  for i = 1:1:length(inds) topD(:,:,inds(i)) = top.Ds(:,:,2); end
        inds = []; inds = find(Vrs==2);  for i = 1:1:length(inds) topD(:,:,inds(i)) = top.Ds(:,:,3); end
        inds = []; inds = find(Vrs==3);  for i = 1:1:length(inds) topD(:,:,inds(i)) = top.Ds(:,:,4); end
        inds = []; inds = find(Vrs==-1); for i = 1:1:length(inds) topD(:,:,inds(i)) = top.Ds(:,:,5); end
        inds = []; inds = find(Vrs==-2); for i = 1:1:length(inds) topD(:,:,inds(i)) = top.Ds(:,:,6); end
        inds = []; inds = find(Vrs==-3); for i = 1:1:length(inds) topD(:,:,inds(i)) = top.Ds(:,:,7); end
end
% --------------------------------------------------------------


% ---------------------- Build Waveform ------------------------
if(1)
    pers = 2; sam = 1000;                       % periods to show. samples/period.
    time = linspace(0,Ts*pers,sam*pers);        % time vector to plot 
    prevJ = 1; j = 1; x0 = Xss; diffIs =[];     % initializations
    for i = 1:1:length(time)
        tI = mod(time(i),Ts); 
        j = find(tI<te,1);
        if(isempty(j)) j = length(te); end
        if(j~=1) tI = tI-te(j-1); end
        if(prevJ ~= j && i ~= 1) 
            x0 = wx(:,i-1);  diffIs = [diffIs, i]; end
        wx(:,i) = expm(topA(:,:,j)*tI)*Xss_all(:,j)+topA(:,:,j)\...
            (expm(topA(:,:,j)*tI)-eye(ns))*topB(:,:,j)*U;  
        wVrec(i) = topC(1,:,j)*wx(:,i) + topD(1,:,j)*U;
        for kk = 1:meters
            wmets(i,kk) = topC(kk,:,j)*wx(:,i); end
        wVinv(i) = Us(j)*U(1); 
        prevJ = j; 
    end    
    
%     figure(2); hold off; 
%     plot(time*1e6, wVrec, 'lineWidth', 2); hold on; 
%     [~, diffIs] = ismember(time(abs(diff(wVrec))>2),time); 
%         diffIs = diffIs(4:15); 
%         diffIs = repelem(diffIs,2); 
%         diffIs(2:2:end) = diffIs(2:2:end)+1;
%         diffVrecs = abs(diff(wVrec(diffIs))); 
%         diffVrecs(diffVrecs<2) = []; 
%     plot(time(diffIs)*1e6, wVrec(diffIs), 'or', 'lineWidth', 2); 
%     hold off; plot(diffVrecs, 'lineWidth', 2); 
% %     plot(time*1e6, meter2, 'lineWidth', 2); 
% %     plot(time*1e6, meter3, 'lineWidth', 2); 
    
%     figure(3); hold off; 
%     plot(time*1e6, wVrec, 'lineWidth', 2); hold on; 
%     plot(time(diffIs)*1e6, wVrec(diffIs), 'or', 'lineWidth', 2); 
%     plot(time*1e6, -wx(indIs,:)*5, 'lineWidth', 2);    
%     xlim([0, 2*Ts*1e6]); 
    
%     ff = figure(2); hold off; font = 'Times'; fntSize = 15;
%     set(ff, 'DefaultTextFontName', font, 'DefaultAxesFontName', font, 'position', [130+sc2 570 860 420]); 
%     if(meters>=6)
%         subplot(3,1,1); hold off; 
%             plot(time*1e6, wmets(:,2), 'lineWidth', 2); hold on; 
%             plot(time*1e6, wmets(:,6), 'lineWidth', 2); 
%             ylabel('Charge Share DS Voltage'); 
%             legend('Low [meter: 2]', 'High [meter: 6]'); end
%     subplot(3,1,2); hold off;
%         plot(time*1e6, wx(2,:), 'lineWidth', 2); hold on; 
%         plot(time*1e6, wx(3,:), 'lineWidth', 2); 
%         plot(time*1e6, wx(indVout,:), 'lineWidth', 4); 
%         ylabel('Flying Cap Voltages'); 
% %         legend('CA1 (High) [state: 2]', 'CA2 (Mid) [state: 3]', 'CA3 (GND) [state: 4]', 'Cout'); 
%         legend('C_{A1}', 'C_{A2}', 'C_{A3}, C_{out}'); 
%     if(meters >= 5)
%         subplot(3,1,3); hold off; 
%             plot(time*1e6, wmets(:,3), 'lineWidth', 2); hold on; 
%             plot(time*1e6, wmets(:,4), 'lineWidth', 2); 
%             plot(time*1e6, wmets(:,5), 'lineWidth', 2); 
%             ylabel('Amps'); 
%             legend('Low [meter: 3]', 'SW node [meter: 4]', 'High [meter: 5]'); end
        
    ff = figure(2); hold off; font = 'Times'; fntSize = 15;
    set(ff, 'DefaultTextFontName', font, 'DefaultAxesFontName', font, 'position', [130+sc2 570 760 350]); 
        hh1 = plot(time*1e6, wx(2,:), 'lineWidth', 2); hold on; 
        hh2 = plot(time*1e6, wx(3,:), 'lineWidth', 2); 
        hh3 = plot(time*1e6, wx(indVout,:), 'lineWidth', 4); 
        ylim([0.95*min([wx(2,:),wx(3,:),wx(indVout)]),1.05*max([wx(2,:),wx(3,:),wx(indVout)])]);
        xlim([0,time(end)*1e6]); 
        ylabel('Voltage, [V]'); 
        xlabel('Time [\muS]'); 
            ax = ancestor(hh1, 'axes'); 
            xrule = ax.XAxis; xrule.FontSize = fntSize; 
            yrule = ax.YAxis(1); yrule.FontSize = fntSize;    
            ll = legend([hh1,hh2,hh3], '{\it{v}}_{{\it{CA}}1}', '{\it{v}}_{{\it{CA}}2}', '{\it{v}}_{{\it{CA}}3}, {\it{v}}_{{\it{Cout}}}'); 
            ll.FontSize = fntSize-1;
            ll.FontSize = fntSize-1;

%     figure(3); hold off; 
% %         plot(time*1e6, wmets(:,3), 'lineWidth', 2); hold on; 
%         plot(time*1e6, wmets(:,4), 'lineWidth', 2); hold on; 
%         plot(time*1e6, wmets(:,5), 'lineWidth', 2); 
%         plot(time*1e6, wmets(:,7), 'lineWidth', 2); 
%         plot(time*1e6, wmets(:,8), 'lineWidth', 2); 
%         plot(time*1e6, wmets(:,9), 'lineWidth', 2);  
%         plot(time*1e6, wmets(:,10), 'lineWidth', 2); 
%         plot(time*1e6, wmets(:,11), 'lineWidth', 2); 
%         ylabel('Amps'); 
%         legend('CS: SW', 'CS: High',...
%                 'CB3', 'CoutN', 'CoutF', 'CoutN1', 'CA3'); 
        
    
%     ff = figure(2); hold off; font = 'Times'; fntSize = 15;
%     set(ff, 'DefaultTextFontName', font, 'DefaultAxesFontName', font, 'position', [130+sc2 570 860 420]);
%         hh1 = plot(time*1e6, wx(indVout,:), 'lineWidth', 2); 
%         ylim([4.99, 5.01]); 
%         xlim([0, max(time*1e6)]); 
%         ylabel('Volage  [V]'); 
%             ax = ancestor(hh1, 'axes'); 
%             xrule = ax.XAxis; xrule.FontSize = fntSize; 
%             yrule = ax.YAxis(1); yrule.FontSize = fntSize;    
%             ll = legend([hh1], '{\it v_{Cout}}'); 
%             ll.FontSize = fntSize-0;
%             ll.FontSize = fntSize-0;


    ff = figure(1); hold off; font = 'Times'; fntSize = 15;
    set(ff, 'DefaultTextFontName', font, 'DefaultAxesFontName', font, 'position', [130+sc2 570 860 420]);
        subplot(2,1,1); hold off; 
        hh1 = plot(time*1e6, wVinv, 'lineWidth', 2); hold on;  
        hh2 = plot(time*1e6, wVrec, 'lineWidth', 2); 
        hh3 = plot(time*1e6, wx(indCp,:), 'lineWidth', 2);
        hh4 = plot(time*1e6, wx(indCs,:), 'lineWidth', 2);     
        ylim([-abs(min(wx(indCs,:)))-0.1.*abs(min(wx(indCs,:))), 1.1*max(wx(indCs,:))]); 
        xlim([0, max(time*1e6)]); 
        ylabel('Volage  [V]'); 
            ax = ancestor(hh1, 'axes'); 
            xrule = ax.XAxis; xrule.FontSize = fntSize; 
            yrule = ax.YAxis(1); yrule.FontSize = fntSize;    
            ll = legend([hh1,hh2,hh3,hh4], '{\it v_{inv}}', '{\it v_{rec}}', '{\it v_{Cp}}', '{\it v_{Cs}}');
            ll.FontSize = fntSize-0;
            ll.FontSize = fntSize-0;
    subplot(2,1,2); hold off; 
        hh1 = plot(time*1e6, wx(indIp,:), 'lineWidth', 2); hold on; 
        hh2 = plot(time*1e6, -wx(indIs,:), 'lineWidth', 2);           % PLECS polarity wrong.
        xlabel('Time [{\mu}s]'); ylabel('Current  [A]'); 
        ylim([-abs(min(wx(indIs,:)))-0.1.*abs(min(wx(indIs,:))), 1.1*max(wx(indIs,:))]); 
        xlim([0, max(time*1e6)]); 
            ax = ancestor(hh1, 'axes'); 
            xrule = ax.XAxis; xrule.FontSize = fntSize; 
            yrule = ax.YAxis(1); yrule.FontSize = fntSize;    
            ll = legend([hh1,hh2], '{\it i_{inv}}', '{\it i_{rec}}'); 
            ll.FontSize = fntSize-0;
            ll.FontSize = fntSize-0;
%     figure(2); hold off; 
%         plot(time*1e6, wx(2,:), 'lineWidth', 2); hold on; 
%         plot(time*1e6, wx(3,:), 'lineWidth', 2);
%         plot(time*1e6, wx(4,:), 'lineWidth', 2);              
%         legend('CA1', 'CA2', 'CA3');
%         ylabel('Volage  [V]'); 
end
% --------------------------------------------------------------


% FROM Rload CODE -- will ensure all is correct before deleting.
% -------------------- Vcs Zero Crossing ----------------------- 
if(0)
    xtn = @(t,n) expm(topA(:,:,n)*t)*Xss_all(:,n)+ ...                  % States in interval 'n' from 
        topA(:,:,n)\(expm(topA(:,:,n)*t)-eye(ns))*topB(:,:,n)*U;     	% ... 'Xss' to 't'.
    vcsISO = 0.*top.Bs(:,:,1)'; vcsISO(indCs) = 1; 
    Vcstn = @(t,n) vcsISO * xtn(t,n);                                   % Vcs only. 
    for zint = 1:1:length(topA(1,1,:))                                  % Find interval where Vcsz occurs.
        if(sign(Vcstn(td(zint),zint)) ~= sign(Xss(indCs))) break; end; end
    LB = 0; UB = td(zint);         
    Vcst = @(t) Vcstn(t,zint);                                          % Vcs(t) w/in 'n' for lsqnonlin.
    options = optimoptions('lsqnonlin', 'Display', 'off');              % Supress display on lsqnonlin.
    tz = lsqnonlin(Vcst,LB,LB,UB,options);                              % Time of zero crossing. minimizing: vcs_II_t^2.
    Xz = xtn(tz,zint);                                                  % States at 'tz' w/in interval 'n'.
    Vcsz_d = vcsISO * (topA(:,:,zint)*Xz + topB(:,:,zint)*U);           % Slope at Vcs(tz,n).  
    phRef = round(-(sum(td(1:zint-1))+tz)/Pclk150); 

    if(0) % USE THIS FIGURE TO CHECK IF THE ZERO CROSSING IS CORRECT!! 
    figure(2); hold off; 
    plot(time*1e6, wVinv, 'lineWidth', 2); hold on;  
    plot(time*1e6, wVrec, 'lineWidth', 2); 
    plot(time*1e6, wx(indCp,:), 'lineWidth', 2);
    plot(time*1e6, wx(indCs,:), 'lineWidth', 2);     
    ylabel('Volage  [V]'); 
    xlim([0 time(end/2)*1e6]); 
    plot((te2)*1e6,Xss_all(indCs,:),'ok'); 
    plot((te2(zint)+tz)*1e6,Xz(indCs),'or','lineWidth',3); 
    plot((te2(zint)+LB)*1e6,Vcst(LB),'ob','lineWidth',3); 
    plot((te2(zint)+UB)*1e6,Vcst(UB),'ob','lineWidth',3);               
    legend('Inv.', 'Rec.', 'Cp', 'Cs', 'Xss', 'Xz', 'zint_0', 'zint_{end}'); end
end
% -------------------------------------------------------------- 


% -------------------- Vcs Zero Crossing ----------------------- 
if(1)
    xtn = @(t,n) expm(topA(:,:,n)*t)*Xss_all(:,n)+ ...                  % States in interval 'n' from 
        topA(:,:,n)\(expm(topA(:,:,n)*t)-eye(ns))*topB(:,:,n)*U;     	% ... 'Xss' to 't'.
    vcsISO = 0.*top.Bs(:,1,1)'; vcsISO(indCs) = 1; 
    Vcstn = @(t,n) vcsISO * xtn(t,n);                                   % Vcs only. 
    for zint = 1:1:length(topA(1,1,:))                                  % Find interval where Vcsz occurs.
        if(sign(Vcstn(td(zint),zint)) ~= sign(Xss(indCs))) break; end; end
    LB = 0; UB = td(zint);         
    Vcst = @(t) Vcstn(t,zint);                                          % Vcs(t) w/in 'n' for lsqnonlin.
    options = optimoptions('lsqnonlin', 'Display', 'off');              % Supress display on lsqnonlin.
    tz = lsqnonlin(Vcst,LB,LB,UB,options);                              % Time of zero crossing. minimizing: vcs_II_t^2.
    Xz = xtn(tz,zint);                                                  % States at 'tz' w/in interval 'n'.
    Vcsz_d = vcsISO * (topA(:,:,zint)*Xz + topB(:,:,zint)*U);           % Slope at Vcs(tz,n).  
    phRef = round(-(sum(td(1:zint-1))+tz)/Pclk150); 

    if(0) % USE THIS FIGURE TO CHECK IF THE ZERO CROSSING IS CORRECT!! 
    figure(2); hold off; 
    plot(time*1e6, wVinv, 'lineWidth', 2); hold on;  
    plot(time*1e6, wVrec, 'lineWidth', 2); 
    plot(time*1e6, wx(indCp,:), 'lineWidth', 2);
    plot(time*1e6, wx(indCs,:), 'lineWidth', 2);     
    ylabel('Volage  [V]'); 
    xlim([0 time(end/2)*1e6]); 
    plot((te2)*1e6,Xss_all(indCs,:),'ok'); 
    plot((te2(zint)+tz)*1e6,Xz(indCs),'or','lineWidth',3); 
    plot((te2(zint)+LB)*1e6,Vcst(LB),'ob','lineWidth',3); 
    plot((te2(zint)+UB)*1e6,Vcst(UB),'ob','lineWidth',3);               
    legend('Inv.', 'Rec.', 'Cp', 'Cs', 'Xss', 'Xz', 'zint_0', 'zint_{end}'); end
end
% -------------------------------------------------------------- 


% ---------------------- Switching Loss ------------------------ 
if(1)
    swCounter = 0; 
    Psw = zeros(size(1:8));      
    for swA = 2:1:8                	% 1/2 peroid of REC actions + 1 INV action.
        if(Us(swA)~=Us(swA-1))   	% INVERTER
            tdc = Pclk150*2; 
            IiT = Coss*Vinv*2/tdc;                                                                       
            Ii = Xss_all(indIp,swA);
            Iia = abs(Ii); 
            if(Ii<0)                                                  
                PswFI = 'B';    
                PswI = Vinv.*f.*(Iia.*trr+Qrr) + Iia.*Vd.*tdc.*f;     
            elseif(Ii>IiT)                                       
                PswFI = 'L';    
                PswI = Iia.*Vd.*f.*(tdc.*(Iia-IiT)./Iia);             
            else                                                  
                PswFI = 'S';
                PswI = 1/2.*Iia.*ton.*f.*(Vinv.*(1-Iia./IiT)) + ...  	 
                     1/2.*Coss.*f.*(Vinv.*(1-Iia./IiT)); end
            PswI = (PswI + Vinv.*Ii.*toff.*f) .* 2;              	
        else                            % RECTIFER 
            swCounter = swCounter + 1; 
            tdc = Pclk150;              % Dead time (diode conduction nominal).
            IiT = Coss*Vout/tdc;       	% Current threshold.              
            Ii = -Xss_all(indIs,swA);  	% PLECS polarity inverted.                                
            Iia = abs(Ii); 
            if(swCounter < 4)                                           % RISING EDGE
                if(Ii<0)                                                % If backwards current.
                    PswF(swA) = 'B';    
                    SW = Vout.*f.*(Iia.*trr+Qrr) + Iia.*Vd.*tdc.*f;     % RR + Diode Cond. 
                elseif(Ii>IiT)                                          % If large current.
                    PswF(swA) = 'L';    
                    SW = Iia.*Vd.*f.*(tdc.*(Iia-IiT)./Iia);             % Diode Cond. 
                else                                                    % If small current.
                    PswF(swA) = 'S';
                    SW = 1/2.*Iia.*ton.*f.*(Vout.*(1-Iia./IiT)) + ...  	% Coss + Ton. 
                         1/2.*Coss.*f.*(Vout.*(1-Iia./IiT)); end
                SW = (SW + Vout.*Ii.*toff.*f) .* 2;                     % Turn off loss, & 1/2 cycle sim.
                Psw(swA) = SW; 
            else
                if(Ii>0)        
                    PswF(swA) = 'B';    
                    SW = Vout.*f.*(Iia.*trr+Qrr) + Iia.*Vd.*tdc.*f;  
                elseif(Ii<-IiT)  
                    PswF(swA) = 'L';    
                    SW = Iia.*Vd.*f.*(tdc.*(Iia-IiT)./Iia);          
                else
                    PswF(swA) = 'S';
                    SW = 1/2.*Iia.*ton.*f.*(Vout.*(1-Iia./IiT)) + ...  	
                         1/2.*Coss.*f.*(Vout.*(1-Iia./IiT)); end
                SW = (SW + Vout.*Ii.*toff.*f) .* 2;                              
                Psw(swA) = SW; 
            end
        end

    end
    lvls = 3; 
        if(M3>=249) Psw(3:4) = 0; end
        if(M2>=248) Psw(2:5) = 0; end
        if(M1>=247) Psw(1:6) = 0; end
end
% -------------------------------------------------------------- 

% ------------------------ Fundamental ------------------------- 
if(1)
    timeLoc = linspace(0,Ts*pers*100,sam*pers*100);
    vrecLoc = repmat(wVrec,1,100); 
    vinvLoc = repmat(wVinv,1,100); 
    irecLoc = repmat(-wx(indIs,:),1,100); 
    iinvLoc = repmat(wx(indIp,:),1,100); 

    [VsVpPhi, vrecFund, vinvFund, THDvrec, THDvinv] = wavesPhi(timeLoc,vrecLoc,vinvLoc,150e3,0); 
    [recPhi,  vrecFund, irecFund, THDvrec, THDirec] = wavesPhi(timeLoc,vrecLoc,irecLoc,150e3,0); 
    [invPhi,  vinvFund, iinvFund, THDvinv, THDiinv] = wavesPhi(timeLoc,vinvLoc,iinvLoc,150e3,0); 
    ZrecFund = vrecFund/irecFund; 
    ZrecPhi = recPhi; 


    [recPhi,  ~,~,~,~] = wavesPhi(timeLoc, irecLoc, vrecLoc, 150e3,0) 
end
% -------------------------------------------------------------- 
 
% ----------- Evaluate Power + Print Operating Point -----------                             
if(1)
    Pin = sum(Us.*U(1).*Xact(indIp,:).*td/Ts);
    if(~exist('Rload')) Rload = Vout/Iload; end
    Pout = Vout^2/Rload - PswI - sum(Psw);               % How I've been calc'ing it. 
    eff = Pout/Pin * 100; 
    
    Ohm = char(937);   
    fprintf( 'Vin:                 %7.2f  [V]\n', Vinv); 
    fprintf( 'Vout:                %7.2f  [V]\n', Vout); 
    fprintf(['Rout:                %7.2f  [' Ohm ']\n'], Rload); 
    fprintf( 'Iout:                %7.2f  [A]\n', Vout/Rload);  
    fprintf( 'Pin:                 %7.2f  [W]\n', Pin); 
    fprintf( 'Pout:                %7.2f  [W]\n', Pout);
    fprintf( 'Eff:                 %7.2f  [%%]\n', eff); 
    fprintf( 'Phase [Vrec|Vcsz]:   %7.0f  [clk]\n', phRef); 
    fprintf( 'Phase [Vinv|Vrec]:   %7.0f  [clk]\n', ph); 
    fprintf( 'Duty:                %7.0f  [%%]\n', round(duty)); 
    fprintf( 'Modulation:          %7.0f  [bit]\n', b11Mod); 
    fprintf( 'M1:                  %7.0f  [clk]\n', M1); 
    fprintf( 'M2:                  %7.0f  [clk]\n', M2); 
    fprintf( 'M3:                  %7.0f  [clk]\n', M3); 
    disp(['------------------------------------']);     
    fprintf( 'Vrec,1 Mag:          %7.2f  [V]\n', vrecFund); 
    fprintf( 'Irec,1 Mag:          %7.2f  [A]\n', irecFund); 
    fprintf( 'Phase [Vinv|Vrec]:   %7.2f  [deg]\n', VsVpPhi); 
    fprintf(['Zrec,1 Mag:          %7.2f  [' Ohm ']\n'], ZrecFund); 
    fprintf( 'Phase [Zrec,1]:      %7.2f  [deg]\n', ZrecPhi); 
    disp(['------------------------------------']); 
    fprintf( 'VCA1:                %7.4f  [V]\n', Vca1); 
    fprintf( 'VCA2:                %7.4f  [V]\n', Vca2);  
    fprintf( 'VCA3:                %7.4f  [V]\n', Vca3);  
    fprintf( 'VCB3:                %7.4f  [V]\n', Vcb3); 
    fprintf( 'VCB2:                %7.4f  [V]\n', Vcb2);
    fprintf( 'VCB1:                %7.4f  [V]\n', Vcb1); 
    disp(['------------------------------------']);   
end
% -------------------------------------------------------------- 

        
%% --------------- Individual Level Pert Info ------------------ 
if(1)
    % NOTE: sqeuence ->   M1    M2    M3
    %                    SA3   SA2   SA1
    %              LVL:   3     2     1
    %                    Low   Mid   High
    inPert = [3,2,1,1,2,3]; inPert = [inPert, inPert]';             % Defines the LVL of each pert.
    lvl = 3; % (LOWEST)
        row1 = (inPert==lvl) .* inPair(:,1); row1(row1==0) = []; 
        row2 = (inPert==lvl) .* inPair(:,2); row2(row2==0) = []; 
        inPair3 = [row1 row2]; 
        inPolD3 = inPolD(inPert==lvl); 
        inPolP3 = inPolP(inPert==lvl); 
    lvl = 2; % (MIDDLE)
        row1 = (inPert==lvl) .* inPair(:,1); row1(row1==0) = []; 
        row2 = (inPert==lvl) .* inPair(:,2); row2(row2==0) = []; 
        inPair2 = [row1 row2]; 
        inPolD2 = inPolD(inPert==lvl); 
        inPolP2 = inPolP(inPert==lvl); 
    lvl = 1; % (HIGHEST)
        row1 = (inPert==lvl) .* inPair(:,1); row1(row1==0) = []; 
        row2 = (inPert==lvl) .* inPair(:,2); row2(row2==0) = []; 
        inPair1 = [row1 row2]; 
        inPolD1 = inPolD(inPert==lvl); 
        inPolP1 = inPolP(inPert==lvl); 

    % Changes overall perterbation info if any LVL is zeroed out. 
    if(M2==250)                             % Remove lvls 1 + 2
        inPair = inPair3; 
        inPolD = inPolD3; 
        inPolP = inPolP3; 
    elseif(M3==250)                         % Remove lvl 1
        inPair(inPert==1,:) = 0;      
        inPair(inPair(:,1)==0,:) = []; 
        inPolD(inPert==1) = []; 
        inPolP(inPert==1) = [];   
    end
        
end
% --------------------------------------------------------------  


% FROM Rload CODE -- will ensure all is correct before deleting. 
% ----------- Transfer Functions (Full Cycle Model) ------------
% ------------------------ (All Levels) ------------------------
if(0)
    % Natural Response Matrix.
    N = 1;                                               
    for jj = fliplr(1:1:length(td))                     % From last interval to 1st.
        N = N * expm(topA(:,:,jj)*td(jj)); end          % propogate through "reverse" period.

    % Forced Response DUTY Matrix.
    Fd = zeros(size([Xss Xss]));                      	
    for ii = 1:1:length(inPair(:,1))
        in1 = inPair(ii,1); in2 = inPair(ii,2);         % intervals b4+after perterbation.
        Xp = Xss_all(:,in2);                            % Xp is steady state op point. 
        if(inPolD(ii) == 1)                            	% "1" is perturbation left to right.
            Ae = topA(:,:,in1); Be = topB(:,in1); 
            As = topA(:,:,in2); Bs = topB(:,in2); 
        else                                            % "0" is perturbation right to left.
            As = topA(:,:,in1); Bs = topB(:,in1); 
            Ae = topA(:,:,in2); Be = topB(:,in2);     
        end
        xd = ((eye(ns)+Ae)*Xp+Be*U) - ...            	% state diff due to perterbation.
                ((eye(ns)+As)*Xp+Bs*U); 
        Fdi = 1;                                      	% initialize the multiplication. 
        for jj = fliplr(in2:1:length(td)) 
            Fdi = Fdi * expm(topA(:,:,jj)*td(jj)); end  
        Fd(:,ii) = Fdi*xd;                             	% all propogated perturbations.
    end
    Fd = sum(Fd,2);                                  	% sum all propogated perturbations.

    % Forced Response PHASE Matrix.
    Fp = zeros(size([Xss Xss]));                      	
    for ii = 1:1:length(inPair(:,1))
        in1 = inPair(ii,1); in2 = inPair(ii,2);         % intervals b4+after perterbation.
        Xp = Xss_all(:,in2);                            % Xp is steady state op point. 
        if(inPolP(ii) == 1)                            	% "1" is perturbation left to right.
            Ae = topA(:,:,in1); Be = topB(:,in1); 
            As = topA(:,:,in2); Bs = topB(:,in2); 
        else                                            % "0" is perturbation right to left.
            As = topA(:,:,in1); Bs = topB(:,in1); 
            Ae = topA(:,:,in2); Be = topB(:,in2);     
        end
        xd = ((eye(ns)+Ae)*Xp+Be*U) - ...             	% state diff due to perterbation.
                ((eye(ns)+As)*Xp+Bs*U); 
        Fdi = 1;                                      	% initialize the multiplication. 
        for jj = fliplr(in2:1:length(td)) 
            Fdi = Fdi * expm(topA(:,:,jj)*td(jj)); end  
        Fp(:,ii) = Fdi*xd;                             	% all propogated perturbations.
    end
    Fp = sum(Fp,2);                                  	% sum all propogated perturbations.

    % Output: Natural Response Matrix.
    H = expm(topA(:,:,zint)*tz);                      	% Zero-crossing interval eval'ed first. 
    for jj = fliplr(1:1:zint)   
        if(jj~=zint) H = H*expm(topA(:,:,jj)*td(jj)); 	% Back to start of period. 
        end; end                                       

    % Output: Forced Response PHASE Matrix.
    Jp = zeros(size([Xss Xss]));       
    rPIs = (inPair(:,1)<zint); rPIs(rPIs==0) = [];      % Relevant perturbations that  
    rPair = [inPair(rPIs,1) inPair(rPIs,2)];            % ... affect tz in FC Model.
    for ii = 1:1:length(rPair(:,1))
        in1 = rPair(ii,1); in2 = rPair(ii,2);           % intervals b4+after perterbation.
        Xp = Xss_all(:,in2);                            % Xp is steady state op point. 
        if(inPolP(ii) == 1)                            	% "1" is perturbation left to right.
            Ae = topA(:,:,in1); Be = topB(:,in1); 
            As = topA(:,:,in2); Bs = topB(:,in2); 
        else                                            % "0" is perturbation right to left.
            As = topA(:,:,in1); Bs = topB(:,in1); 
            Ae = topA(:,:,in2); Be = topB(:,in2);     
        end
        xd = ((eye(ns)+Ae)*Xp+Be*U) - ...            	% state diff due to perterbation.
                ((eye(ns)+As)*Xp+Bs*U); 
        Fdi = expm(topA(:,:,zint)*tz);               	% initialize the multiplication. 
        for jj = fliplr(in2:1:zint)                     % From pert to zero-cross interval.
            if(jj~=zint) Fdi = Fdi*...
                expm(topA(:,:,jj)*td(jj)); 
            end; end; Jp(:,ii) = Fdi*xd;             	% all propogated perturbations.
    end
    Jp = sum(Jp,2);                                

    % Output: Forced Response PHASE Matrix.
    Jd = zeros(size([Xss Xss]));       
    rPIs = (inPair(:,1)<zint); rPIs(rPIs==0) = [];      % Relevant perturbations that  
    rPair = [inPair(rPIs,1) inPair(rPIs,2)];            % ... affect tz in FC Model.
    for ii = 1:1:length(rPair(:,1))
        in1 = rPair(ii,1); in2 = rPair(ii,2);           % intervals b4+after perterbation.
        Xp = Xss_all(:,in2);                            % Xp is steady state op point. 
        if(inPolD(ii) == 1)                            	% "1" is perturbation left to right.
            Ae = topA(:,:,in1); Be = topB(:,in1); 
            As = topA(:,:,in2); Bs = topB(:,in2); 
        else                                            % "0" is perturbation right to left.
            As = topA(:,:,in1); Bs = topB(:,in1); 
            Ae = topA(:,:,in2); Be = topB(:,in2);     
        end
        xd = ((eye(ns)+Ae)*Xp+Be*U) - ...            	% state diff due to perterbation.
                ((eye(ns)+As)*Xp+Bs*U); 
        Fdi = expm(topA(:,:,zint)*tz);               	% initialize the multiplication. 
        for jj = fliplr(in2:1:zint)                     % From pert to zero-cross interval.
            if(jj~=zint) Fdi = Fdi*...
                expm(topA(:,:,jj)*td(jj)); 
            end; end; Jd(:,ii) = Fdi*xd;             	% all propogated perturbations.
    end
    Jd = sum(Jd,2);                                  	% sum all propogated perturbations.

    if(M1==250) 
        N = N*0; Fd = Fd*0; Fp = Fp*0; 
        H = H*0; Jd = Jd*0; Jp = Jp*0; end
    
    voutISO = 0.*top.Bs(:,:,1)';      	% Isolate state 5 (Vout).
        voutISO(indVout) = 1; 
    Gvd = ss(N,Fd,voutISO,0,Ts);        % Time[duty]  -> Vout: ISO * ((zI-N)^(-1)*Fd)
    Gvp = ss(N,Fp,voutISO,0,Ts);        % Time[phase] -> Vout: ISO * ((zI-N)^(-1)*Fp)
    Gpp = ss(N,Fp,H,Jp,Ts);             % Here,                  H * ((zI-N)^(-1)*Fp) + Jp
        Gpp = Gpp(indCs,1)/(-Vcsz_d);   % Time[phase] -> Time[zero-x]
        Gpp = -1 * Gpp;
    Gpd = ss(N,Fd,H,Jd,Ts);             % Here,                  H * ((zI-N)^(-1)*Fd) + Jd
        Gpd = Gpd(indCs,1)/(-Vcsz_d);   % Time[phase] -> Time[zero-x]

    [GvdN,GvdD] = tfdata(Gvd,'v');  	% Numerator + denominator data of Gvd.
    [GvpN,GvpD] = tfdata(Gvp,'v');  	% Numerator + denominator data of Gvp.
    [GppN,GppD] = tfdata(Gpp,'v');  	% Numerator + denominator data of Gpp.
    [GpdN,GpdD] = tfdata(Gpd,'v');  	% Numerator + denominator data of Gpd.
end
% -------------------------------------------------------------- 


% FROM Rload CODE -- will ensure all is correct before deleting. 
% ----------- Transfer Functions (Full Cycle Model) ------------
% -------------------- (Individual LVLs) -----------------------
if(0)
    % Forced Response DUTY Matrix. (LVL 1)
    Fd1 = zeros(size([Xss Xss]));                      	
    for ii = 1:1:length(inPair1(:,1))
        in1 = inPair1(ii,1); in2 = inPair1(ii,2);   	% intervals b4+after perterbation.
        Xp = Xss_all(:,in2);                            % Xp is steady state op point. 
        if(inPolD1(ii) == 1)                           	% "1" is perturbation left to right.
            Ae = topA(:,:,in1); Be = topB(:,in1); 
            As = topA(:,:,in2); Bs = topB(:,in2); 
        else                                            % "0" is perturbation right to left.
            As = topA(:,:,in1); Bs = topB(:,in1); 
            Ae = topA(:,:,in2); Be = topB(:,in2);     
        end
        xd = ((eye(ns)+Ae)*Xp+Be*U) - ...            	% state diff due to perterbation.
                ((eye(ns)+As)*Xp+Bs*U); 
        Fdi1 = 1;                                      	% initialize the multiplication. 
        for jj = fliplr(in2:1:length(td)) 
            Fdi1 = Fdi1 * expm(topA(:,:,jj)*td(jj)); end  
        Fd1(:,ii) = Fdi1*xd;                           	% all propogated perturbations.
    end
    Fd1 = sum(Fd1,2);                                  	% sum all propogated perturbations.

    % Forced Response DUTY Matrix. (LVL 2)
    Fd2 = zeros(size([Xss Xss]));                      	
    for ii = 1:1:length(inPair2(:,1))
        in1 = inPair2(ii,1); in2 = inPair2(ii,2);   	% intervals b4+after perterbation.
        Xp = Xss_all(:,in2);                            % Xp is steady state op point. 
        if(inPolD2(ii) == 1)                           	% "1" is perturbation left to right.
            Ae = topA(:,:,in1); Be = topB(:,in1); 
            As = topA(:,:,in2); Bs = topB(:,in2); 
        else                                            % "0" is perturbation right to left.
            As = topA(:,:,in1); Bs = topB(:,in1); 
            Ae = topA(:,:,in2); Be = topB(:,in2);     
        end
        xd = ((eye(ns)+Ae)*Xp+Be*U) - ...            	% state diff due to perterbation.
                ((eye(ns)+As)*Xp+Bs*U); 
        Fdi2 = 1;                                      	% initialize the multiplication. 
        for jj = fliplr(in2:1:length(td)) 
            Fdi2 = Fdi2 * expm(topA(:,:,jj)*td(jj)); end  
        Fd2(:,ii) = Fdi2*xd;                           	% all propogated perturbations.
    end
    Fd2 = sum(Fd2,2);                                  	% sum all propogated perturbations.

    % Forced Response DUTY Matrix. (LVL 3)
    Fd3 = zeros(size([Xss Xss]));                      	
    for ii = 1:1:length(inPair3(:,1))
        in1 = inPair3(ii,1); in2 = inPair3(ii,2);   	% intervals b4+after perterbation.
        Xp = Xss_all(:,in2);                            % Xp is steady state op point. 
        if(inPolD3(ii) == 1)                           	% "1" is perturbation left to right.
            Ae = topA(:,:,in1); Be = topB(:,in1); 
            As = topA(:,:,in2); Bs = topB(:,in2); 
        else                                            % "0" is perturbation right to left.
            As = topA(:,:,in1); Bs = topB(:,in1); 
            Ae = topA(:,:,in2); Be = topB(:,in2);     
        end
        xd = ((eye(ns)+Ae)*Xp+Be*U) - ...            	% state diff due to perterbation.
                ((eye(ns)+As)*Xp+Bs*U); 
        Fdi3 = 1;                                      	% initialize the multiplication. 
        for jj = fliplr(in2:1:length(td)) 
            Fdi3 = Fdi3 * expm(topA(:,:,jj)*td(jj)); end  
        Fd3(:,ii) = Fdi3*xd;                           	% all propogated perturbations.
    end
    Fd3 = sum(Fd3,2);                                  	% sum all propogated perturbations.

    % Output: Forced Response DUTY Matrix. (LVL 1)
    Jd1 = zeros(size([Xss Xss]));       
    rPIs = (inPair1(:,1)<zint); rPIs(rPIs==0) = [];    	% Relevant perturbations that  
    rPair = [inPair1(rPIs,1) inPair1(rPIs,2)];          % ... affect tz in FC Model.
    for ii = 1:1:length(rPair(:,1))
        in1 = rPair(ii,1); in2 = rPair(ii,2);           % intervals b4+after perterbation.
        Xp = Xss_all(:,in2);                            % Xp is steady state op point. 
        if(inPolD1(ii) == 1)                         	% "1" is perturbation left to right.
            Ae = topA(:,:,in1); Be = topB(:,in1); 
            As = topA(:,:,in2); Bs = topB(:,in2); 
        else                                            % "0" is perturbation right to left.
            As = topA(:,:,in1); Bs = topB(:,in1); 
            Ae = topA(:,:,in2); Be = topB(:,in2);     
        end
        xd = ((eye(ns)+Ae)*Xp+Be*U) - ...            	% state diff due to perterbation.
                ((eye(ns)+As)*Xp+Bs*U); 
        Jdi1 = expm(topA(:,:,zint)*tz);               	% initialize the multiplication. 
        for jj = fliplr(in2:1:zint)                     % From pert to zero-cross interval.
            if(jj~=zint) Jdi1 = Jdi1*...
                expm(topA(:,:,jj)*td(jj)); 
            end; end; Jd1(:,ii) = Jdi1*xd;             	% all propogated perturbations.
    end
    Jd1 = sum(Jd1,2);                                  	% sum all propogated perturbations.

    % Output: Forced Response DUTY Matrix. (LVL 2)
    Jd2 = zeros(size([Xss Xss]));       
    rPIs = (inPair2(:,1)<zint); rPIs(rPIs==0) = [];  	% Relevant perturbations that  
    rPair = [inPair2(rPIs,1) inPair2(rPIs,2)];       	% ... affect tz in FC Model.
    for ii = 1:1:length(rPair(:,1))
        in1 = rPair(ii,1); in2 = rPair(ii,2);           % intervals b4+after perterbation.
        Xp = Xss_all(:,in2);                            % Xp is steady state op point. 
        if(inPolD2(ii) == 1)                           	% "1" is perturbation left to right.
            Ae = topA(:,:,in1); Be = topB(:,in1); 
            As = topA(:,:,in2); Bs = topB(:,in2); 
        else                                            % "0" is perturbation right to left.
            As = topA(:,:,in1); Bs = topB(:,in1); 
            Ae = topA(:,:,in2); Be = topB(:,in2);     
        end
        xd = ((eye(ns)+Ae)*Xp+Be*U) - ...            	% state diff due to perterbation.
                ((eye(ns)+As)*Xp+Bs*U); 
        Jdi2 = expm(topA(:,:,zint)*tz);               	% initialize the multiplication. 
        for jj = fliplr(in2:1:zint)                     % From pert to zero-cross interval.
            if(jj~=zint) Jdi2 = Jdi2*...
                expm(topA(:,:,jj)*td(jj)); 
            end; end; Jd2(:,ii) = Jdi2*xd;             	% all propogated perturbations.
    end
    Jd2 = sum(Jd2,2);                                  	% sum all propogated perturbations.

    % Output: Forced Response DUTY Matrix. (LVL 3)
    Jd3 = zeros(size([Xss Xss]));       
    rPIs = (inPair3(:,1)<zint); rPIs(rPIs==0) = [];   	% Relevant perturbations that  
    rPair = [inPair3(rPIs,1) inPair3(rPIs,2)];        	% ... affect tz in FC Model.
    for ii = 1:1:length(rPair(:,1))
        in1 = rPair(ii,1); in2 = rPair(ii,2);           % intervals b4+after perterbation.
        Xp = Xss_all(:,in2);                            % Xp is steady state op point. 
        if(inPolD3(ii) == 1)                          	% "1" is perturbation left to right.
            Ae = topA(:,:,in1); Be = topB(:,in1); 
            As = topA(:,:,in2); Bs = topB(:,in2); 
        else                                            % "0" is perturbation right to left.
            As = topA(:,:,in1); Bs = topB(:,in1); 
            Ae = topA(:,:,in2); Be = topB(:,in2);     
        end
        xd = ((eye(ns)+Ae)*Xp+Be*U) - ...            	% state diff due to perterbation.
                ((eye(ns)+As)*Xp+Bs*U); 
        Jdi3 = expm(topA(:,:,zint)*tz);               	% initialize the multiplication. 
        for jj = fliplr(in2:1:zint)                     % From pert to zero-cross interval.
            if(jj~=zint) Jdi3 = Jdi3*...
                expm(topA(:,:,jj)*td(jj)); 
            end; end; Jd3(:,ii) = Jdi3*xd;             	% all propogated perturbations.
    end
    Jd3 = sum(Jd3,2);                                  	% sum all propogated perturbations.

    Gvd1 = ss(N,Fd1,voutISO,0,Ts);                          % Time[duty]  -> Vout: ISO * ((zI-N)^(-1)*Fd)
        if(M3==250) Gvd1 = ss(N*0,Fd1*0,voutISO,0,Ts); end
    Gvd2 = ss(N,Fd2,voutISO,0,Ts);  
        if(M2==250) Gvd2 = ss(N*0,Fd2*0,voutISO,0,Ts); end
    Gvd3 = ss(N,Fd3,voutISO,0,Ts);  
        if(M1==250) Gvd3 = ss(N*0,Fd3*0,voutISO,0,Ts); end

    Gpd1 = ss(N,Fd1,H,Jd1,Ts);                              % Here, H * ((zI-N)^(-1)*Fd) + Jd
        if(M3==250) Gpd1 = ss(N*0,Fd1*0,H*0,Jd1*0,Ts); end  % Time[phase] -> Time[zero-x]
        Gpd1 = Gpd1(indCs,1)/(-Vcsz_d);                     
    Gpd2 = ss(N,Fd2,H,Jd2,Ts);            
        if(M2==250) Gpd2 = ss(N*0,Fd2*0,H*0,Jd2*0,Ts); end
        Gpd2 = Gpd2(indCs,1)/(-Vcsz_d);  
    Gpd3 = ss(N,Fd3,H,Jd3,Ts);
        if(M1==250) Gpd3 = ss(N*0,Fd3*0,H*0,Jd3*0,Ts); end
        Gpd3 = Gpd3(indCs,1)/(-Vcsz_d);   

    [GvdN1,GvdD1] = tfdata(Gvd1,'v');     	% Numerator + denominator data of Gvd.
        [GvdN2,GvdD2] = tfdata(Gvd2,'v');  	
        [GvdN3,GvdD3] = tfdata(Gvd3,'v');  		
    [GpdN1,GpdD1] = tfdata(Gpd1,'v');   	% Numerator + denominator data of Gpd.
        [GpdN2,GpdD2] = tfdata(Gpd2,'v');  	
        [GpdN3,GpdD3] = tfdata(Gpd3,'v');  	
end
% --------------------------------------------------------------


% ----------- Transfer Functions (Full Cycle Model) ------------
% ------------------------ (All Levels) ------------------------
if(1)
    % Natural Response Matrix.
    N = 1;                                               
    for jj = fliplr(1:1:length(td))                     % From last interval to 1st.
        N = N * expm(topA(:,:,jj)*td(jj)); end          % propogate through "reverse" period.

    % Forced Response DUTY Matrix.
    Fd = zeros(size([Xss Xss]));                      	
    for ii = 1:1:length(inPair(:,1))
        in1 = inPair(ii,1); in2 = inPair(ii,2);         % intervals b4+after perterbation.
        Xp = Xss_all(:,in2);                            % Xp is steady state op point. 
        if(inPolD(ii) == 1)                            	% "1" is perturbation left to right.
            Ae = topA(:,:,in1); Be = topB(:,:,in1); 
            As = topA(:,:,in2); Bs = topB(:,:,in2); 
        else                                            % "0" is perturbation right to left.
            As = topA(:,:,in1); Bs = topB(:,:,in1); 
            Ae = topA(:,:,in2); Be = topB(:,:,in2);     
        end
        xd = ((eye(ns)+Ae)*Xp+Be*U) - ...            	% state diff due to perterbation.
                ((eye(ns)+As)*Xp+Bs*U); 
        Fdi = 1;                                      	% initialize the multiplication. 
        for jj = fliplr(in2:1:length(td)) 
            Fdi = Fdi * expm(topA(:,:,jj)*td(jj)); end  
        Fd(:,ii) = Fdi*xd;                             	% all propogated perturbations.
    end
    Fd = sum(Fd,2);                                  	% sum all propogated perturbations.

    % Forced Response PHASE Matrix.
    Fp = zeros(size([Xss Xss]));                      	
    for ii = 1:1:length(inPair(:,1))
        in1 = inPair(ii,1); in2 = inPair(ii,2);         % intervals b4+after perterbation.
        Xp = Xss_all(:,in2);                            % Xp is steady state op point. 
        if(inPolP(ii) == 1)                            	% "1" is perturbation left to right.
            Ae = topA(:,:,in1); Be = topB(:,:,in1); 
            As = topA(:,:,in2); Bs = topB(:,:,in2); 
        else                                            % "0" is perturbation right to left.
            As = topA(:,:,in1); Bs = topB(:,:,in1); 
            Ae = topA(:,:,in2); Be = topB(:,:,in2);     
        end
        xd = ((eye(ns)+Ae)*Xp+Be*U) - ...             	% state diff due to perterbation.
                ((eye(ns)+As)*Xp+Bs*U); 
        Fdi = 1;                                      	% initialize the multiplication. 
        for jj = fliplr(in2:1:length(td)) 
            Fdi = Fdi * expm(topA(:,:,jj)*td(jj)); end  
        Fp(:,ii) = Fdi*xd;                             	% all propogated perturbations.
    end
    Fp = sum(Fp,2);                                  	% sum all propogated perturbations.

    % Output: Natural Response Matrix.
    H = expm(topA(:,:,zint)*tz);                      	% Zero-crossing interval eval'ed first. 
    for jj = fliplr(1:1:zint)   
        if(jj~=zint) H = H*expm(topA(:,:,jj)*td(jj)); 	% Back to start of period. 
        end; end                                       

    % Output: Forced Response PHASE Matrix.
    Jp = zeros(size([Xss Xss]));       
    rPIs = (inPair(:,1)<zint); rPIs(rPIs==0) = [];      % Relevant perturbations that  
    rPair = [inPair(rPIs,1) inPair(rPIs,2)];            % ... affect tz in FC Model.
    for ii = 1:1:length(rPair(:,1))
        in1 = rPair(ii,1); in2 = rPair(ii,2);           % intervals b4+after perterbation.
        Xp = Xss_all(:,in2);                            % Xp is steady state op point. 
        if(inPolP(ii) == 1)                            	% "1" is perturbation left to right.
            Ae = topA(:,:,in1); Be = topB(:,:,in1); 
            As = topA(:,:,in2); Bs = topB(:,:,in2); 
        else                                            % "0" is perturbation right to left.
            As = topA(:,:,in1); Bs = topB(:,:,in1); 
            Ae = topA(:,:,in2); Be = topB(:,:,in2);     
        end
        xd = ((eye(ns)+Ae)*Xp+Be*U) - ...            	% state diff due to perterbation.
                ((eye(ns)+As)*Xp+Bs*U); 
        Fdi = expm(topA(:,:,zint)*tz);               	% initialize the multiplication. 
        for jj = fliplr(in2:1:zint)                     % From pert to zero-cross interval.
            if(jj~=zint) Fdi = Fdi*...
                expm(topA(:,:,jj)*td(jj)); 
            end; end; Jp(:,ii) = Fdi*xd;             	% all propogated perturbations.
    end
    Jp = sum(Jp,2);                                

    % Output: Forced Response PHASE Matrix.
    Jd = zeros(size([Xss Xss]));       
    rPIs = (inPair(:,1)<zint); rPIs(rPIs==0) = [];      % Relevant perturbations that  
    rPair = [inPair(rPIs,1) inPair(rPIs,2)];            % ... affect tz in FC Model.
    for ii = 1:1:length(rPair(:,1))
        in1 = rPair(ii,1); in2 = rPair(ii,2);           % intervals b4+after perterbation.
        Xp = Xss_all(:,in2);                            % Xp is steady state op point. 
        if(inPolD(ii) == 1)                            	% "1" is perturbation left to right.
            Ae = topA(:,:,in1); Be = topB(:,:,in1); 
            As = topA(:,:,in2); Bs = topB(:,:,in2); 
        else                                            % "0" is perturbation right to left.
            As = topA(:,:,in1); Bs = topB(:,:,in1); 
            Ae = topA(:,:,in2); Be = topB(:,:,in2);     
        end
        xd = ((eye(ns)+Ae)*Xp+Be*U) - ...            	% state diff due to perterbation.
                ((eye(ns)+As)*Xp+Bs*U); 
        Fdi = expm(topA(:,:,zint)*tz);               	% initialize the multiplication. 
        for jj = fliplr(in2:1:zint)                     % From pert to zero-cross interval.
            if(jj~=zint) Fdi = Fdi*...
                expm(topA(:,:,jj)*td(jj)); 
            end; end; Jd(:,ii) = Fdi*xd;             	% all propogated perturbations.
    end
    Jd = sum(Jd,2);                                  	% sum all propogated perturbations.

    if(M1==250) 
        N = N*0; Fd = Fd*0; Fp = Fp*0; 
        H = H*0; Jd = Jd*0; Jp = Jp*0; end
    
    voutISO = 0.*top.Bs(:,:,1)';      	% Isolate state 5 (Vout).
        voutISO(indVout) = 1; 
    Gvd = ss(N,Fd,voutISO,0,Ts);        % Time[duty]  -> Vout: ISO * ((zI-N)^(-1)*Fd)
        Gvd = Gvd(1,1); 
    Gvp = ss(N,Fp,voutISO,0,Ts);        % Time[phase] -> Vout: ISO * ((zI-N)^(-1)*Fp)
        Gvp = Gvp(1,1); 
    Gpp = ss(N,Fp,H,Jp,Ts);             % Here,                  H * ((zI-N)^(-1)*Fp) + Jp
        Gpp = Gpp(indCs,1)/(-Vcsz_d);   % Time[phase] -> Time[zero-x]
        Gpp = -1 * Gpp;
    Gpd = ss(N,Fd,H,Jd,Ts);             % Here,                  H * ((zI-N)^(-1)*Fd) + Jd
        Gpd = Gpd(indCs,1)/(-Vcsz_d);   % Time[phase] -> Time[zero-x]

    [GvdN,GvdD] = tfdata(Gvd,'v');  	% Numerator + denominator data of Gvd.
    [GvpN,GvpD] = tfdata(Gvp,'v');  	% Numerator + denominator data of Gvp.
    [GppN,GppD] = tfdata(Gpp,'v');  	% Numerator + denominator data of Gpp.
    [GpdN,GpdD] = tfdata(Gpd,'v');  	% Numerator + denominator data of Gpd.
end
% -------------------------------------------------------------- 


% ----------- Transfer Functions (Full Cycle Model) ------------
% -------------------- (Individual LVLs) -----------------------
if(1)
    % Forced Response DUTY Matrix. (LVL 1)
    Fd1 = zeros(size([Xss Xss]));                      	
    for ii = 1:1:length(inPair1(:,1))
        in1 = inPair1(ii,1); in2 = inPair1(ii,2);   	% intervals b4+after perterbation.
        Xp = Xss_all(:,in2);                            % Xp is steady state op point. 
        if(inPolD1(ii) == 1)                           	% "1" is perturbation left to right.
            Ae = topA(:,:,in1); Be = topB(:,:,in1); 
            As = topA(:,:,in2); Bs = topB(:,:,in2); 
        else                                            % "0" is perturbation right to left.
            As = topA(:,:,in1); Bs = topB(:,:,in1); 
            Ae = topA(:,:,in2); Be = topB(:,:,in2);     
        end
        xd = ((eye(ns)+Ae)*Xp+Be*U) - ...            	% state diff due to perterbation.
                ((eye(ns)+As)*Xp+Bs*U); 
        Fdi1 = 1;                                      	% initialize the multiplication. 
        for jj = fliplr(in2:1:length(td)) 
            Fdi1 = Fdi1 * expm(topA(:,:,jj)*td(jj)); end  
        Fd1(:,ii) = Fdi1*xd;                           	% all propogated perturbations.
    end
    Fd1 = sum(Fd1,2);                                  	% sum all propogated perturbations.

    % Forced Response DUTY Matrix. (LVL 2)
    Fd2 = zeros(size([Xss Xss]));                      	
    for ii = 1:1:length(inPair2(:,1))
        in1 = inPair2(ii,1); in2 = inPair2(ii,2);   	% intervals b4+after perterbation.
        Xp = Xss_all(:,in2);                            % Xp is steady state op point. 
        if(inPolD2(ii) == 1)                           	% "1" is perturbation left to right.
            Ae = topA(:,:,in1); Be = topB(:,:,in1); 
            As = topA(:,:,in2); Bs = topB(:,:,in2); 
        else                                            % "0" is perturbation right to left.
            As = topA(:,:,in1); Bs = topB(:,:,in1); 
            Ae = topA(:,:,in2); Be = topB(:,:,in2);     
        end
        xd = ((eye(ns)+Ae)*Xp+Be*U) - ...            	% state diff due to perterbation.
                ((eye(ns)+As)*Xp+Bs*U); 
        Fdi2 = 1;                                      	% initialize the multiplication. 
        for jj = fliplr(in2:1:length(td)) 
            Fdi2 = Fdi2 * expm(topA(:,:,jj)*td(jj)); end  
        Fd2(:,ii) = Fdi2*xd;                           	% all propogated perturbations.
    end
    Fd2 = sum(Fd2,2);                                  	% sum all propogated perturbations.

    % Forced Response DUTY Matrix. (LVL 3)
    Fd3 = zeros(size([Xss Xss]));                      	
    for ii = 1:1:length(inPair3(:,1))
        in1 = inPair3(ii,1); in2 = inPair3(ii,2);   	% intervals b4+after perterbation.
        Xp = Xss_all(:,in2);                            % Xp is steady state op point. 
        if(inPolD3(ii) == 1)                           	% "1" is perturbation left to right.
            Ae = topA(:,:,in1); Be = topB(:,:,in1); 
            As = topA(:,:,in2); Bs = topB(:,:,in2); 
        else                                            % "0" is perturbation right to left.
            As = topA(:,:,in1); Bs = topB(:,:,in1); 
            Ae = topA(:,:,in2); Be = topB(:,:,in2);     
        end
        xd = ((eye(ns)+Ae)*Xp+Be*U) - ...            	% state diff due to perterbation.
                ((eye(ns)+As)*Xp+Bs*U); 
        Fdi3 = 1;                                      	% initialize the multiplication. 
        for jj = fliplr(in2:1:length(td)) 
            Fdi3 = Fdi3 * expm(topA(:,:,jj)*td(jj)); end  
        Fd3(:,ii) = Fdi3*xd;                           	% all propogated perturbations.
    end
    Fd3 = sum(Fd3,2);                                  	% sum all propogated perturbations.

    % Output: Forced Response DUTY Matrix. (LVL 1)
    Jd1 = zeros(size([Xss Xss]));       
    rPIs = (inPair1(:,1)<zint); rPIs(rPIs==0) = [];    	% Relevant perturbations that  
    rPair = [inPair1(rPIs,1) inPair1(rPIs,2)];          % ... affect tz in FC Model.
    for ii = 1:1:length(rPair(:,1))
        in1 = rPair(ii,1); in2 = rPair(ii,2);           % intervals b4+after perterbation.
        Xp = Xss_all(:,in2);                            % Xp is steady state op point. 
        if(inPolD1(ii) == 1)                         	% "1" is perturbation left to right.
            Ae = topA(:,:,in1); Be = topB(:,:,in1); 
            As = topA(:,:,in2); Bs = topB(:,:,in2); 
        else                                            % "0" is perturbation right to left.
            As = topA(:,:,in1); Bs = topB(:,:,in1); 
            Ae = topA(:,:,in2); Be = topB(:,:,in2);     
        end
        xd = ((eye(ns)+Ae)*Xp+Be*U) - ...            	% state diff due to perterbation.
                ((eye(ns)+As)*Xp+Bs*U); 
        Jdi1 = expm(topA(:,:,zint)*tz);               	% initialize the multiplication. 
        for jj = fliplr(in2:1:zint)                     % From pert to zero-cross interval.
            if(jj~=zint) Jdi1 = Jdi1*...
                expm(topA(:,:,jj)*td(jj)); 
            end; end; Jd1(:,ii) = Jdi1*xd;             	% all propogated perturbations.
    end
    Jd1 = sum(Jd1,2);                                  	% sum all propogated perturbations.

    % Output: Forced Response DUTY Matrix. (LVL 2)
    Jd2 = zeros(size([Xss Xss]));       
    rPIs = (inPair2(:,1)<zint); rPIs(rPIs==0) = [];  	% Relevant perturbations that  
    rPair = [inPair2(rPIs,1) inPair2(rPIs,2)];       	% ... affect tz in FC Model.
    for ii = 1:1:length(rPair(:,1))
        in1 = rPair(ii,1); in2 = rPair(ii,2);           % intervals b4+after perterbation.
        Xp = Xss_all(:,in2);                            % Xp is steady state op point. 
        if(inPolD2(ii) == 1)                           	% "1" is perturbation left to right.
            Ae = topA(:,:,in1); Be = topB(:,:,in1); 
            As = topA(:,:,in2); Bs = topB(:,:,in2); 
        else                                            % "0" is perturbation right to left.
            As = topA(:,:,in1); Bs = topB(:,:,in1); 
            Ae = topA(:,:,in2); Be = topB(:,:,in2);     
        end
        xd = ((eye(ns)+Ae)*Xp+Be*U) - ...            	% state diff due to perterbation.
                ((eye(ns)+As)*Xp+Bs*U); 
        Jdi2 = expm(topA(:,:,zint)*tz);               	% initialize the multiplication. 
        for jj = fliplr(in2:1:zint)                     % From pert to zero-cross interval.
            if(jj~=zint) Jdi2 = Jdi2*...
                expm(topA(:,:,jj)*td(jj)); 
            end; end; Jd2(:,ii) = Jdi2*xd;             	% all propogated perturbations.
    end
    Jd2 = sum(Jd2,2);                                  	% sum all propogated perturbations.

    % Output: Forced Response DUTY Matrix. (LVL 3)
    Jd3 = zeros(size([Xss Xss]));       
    rPIs = (inPair3(:,1)<zint); rPIs(rPIs==0) = [];   	% Relevant perturbations that  
    rPair = [inPair3(rPIs,1) inPair3(rPIs,2)];        	% ... affect tz in FC Model.
    for ii = 1:1:length(rPair(:,1))
        in1 = rPair(ii,1); in2 = rPair(ii,2);           % intervals b4+after perterbation.
        Xp = Xss_all(:,in2);                            % Xp is steady state op point. 
        if(inPolD3(ii) == 1)                          	% "1" is perturbation left to right.
            Ae = topA(:,:,in1); Be = topB(:,:,in1); 
            As = topA(:,:,in2); Bs = topB(:,:,in2); 
        else                                            % "0" is perturbation right to left.
            As = topA(:,:,in1); Bs = topB(:,:,in1); 
            Ae = topA(:,:,in2); Be = topB(:,:,in2);     
        end
        xd = ((eye(ns)+Ae)*Xp+Be*U) - ...            	% state diff due to perterbation.
                ((eye(ns)+As)*Xp+Bs*U); 
        Jdi3 = expm(topA(:,:,zint)*tz);               	% initialize the multiplication. 
        for jj = fliplr(in2:1:zint)                     % From pert to zero-cross interval.
            if(jj~=zint) Jdi3 = Jdi3*...
                expm(topA(:,:,jj)*td(jj)); 
            end; end; Jd3(:,ii) = Jdi3*xd;             	% all propogated perturbations.
    end
    Jd3 = sum(Jd3,2);                                  	% sum all propogated perturbations.

    Gvd1 = ss(N,Fd1,voutISO,0,Ts);                          % Time[duty]  -> Vout: ISO * ((zI-N)^(-1)*Fd)
        if(M3==250) Gvd1 = ss(N*0,Fd1*0,voutISO,0,Ts); end
        Gvd1 = Gvd1(1,1); 
    Gvd2 = ss(N,Fd2,voutISO,0,Ts);  
        if(M2==250) Gvd2 = ss(N*0,Fd2*0,voutISO,0,Ts); end
        Gvd2 = Gvd2(1,1); 
    Gvd3 = ss(N,Fd3,voutISO,0,Ts);  
        if(M1==250) Gvd3 = ss(N*0,Fd3*0,voutISO,0,Ts); end
        Gvd3 = Gvd3(1,1); 

    Gpd1 = ss(N,Fd1,H,Jd1,Ts);                              % Here, H * ((zI-N)^(-1)*Fd) + Jd
        if(M3==250) Gpd1 = ss(N*0,Fd1*0,H*0,Jd1*0,Ts); end  % Time[phase] -> Time[zero-x]
        Gpd1 = Gpd1(indCs,1)/(-Vcsz_d);                     
    Gpd2 = ss(N,Fd2,H,Jd2,Ts);            
        if(M2==250) Gpd2 = ss(N*0,Fd2*0,H*0,Jd2*0,Ts); end
        Gpd2 = Gpd2(indCs,1)/(-Vcsz_d);  
    Gpd3 = ss(N,Fd3,H,Jd3,Ts);
        if(M1==250) Gpd3 = ss(N*0,Fd3*0,H*0,Jd3*0,Ts); end
        Gpd3 = Gpd3(indCs,1)/(-Vcsz_d);   

    [GvdN1,GvdD1] = tfdata(Gvd1,'v');     	% Numerator + denominator data of Gvd.
        [GvdN2,GvdD2] = tfdata(Gvd2,'v');  	
        [GvdN3,GvdD3] = tfdata(Gvd3,'v');  		
    [GpdN1,GpdD1] = tfdata(Gpd1,'v');   	% Numerator + denominator data of Gpd.
        [GpdN2,GpdD2] = tfdata(Gpd2,'v');  	
        [GpdN3,GpdD3] = tfdata(Gpd3,'v');  	
end
% --------------------------------------------------------------


fprintf('Execution Time: %3.0f [min]', round(toc/60));  
fprintf(', %2.0f [s]\n', mod(toc,60));  

% ------------- Plot Small Signal Characteristics -------------- 
if(1)
    if(0)
        ff = figure(2); hold off;                       % Time[DUTY]  -> Vout
        set(ff, 'position', [690+sc2 570 560 420]);     % Time[PHASE] -> Vout
        bode(-Gvd,opts); hold on; bode(Gvp,opts);    
        ylim([-360,0]); legend('Gvd', 'Gvp'); title('All Perturbations'); 

        ff = figure(3); hold off;                       % Time[DUTY]  -> Time[ZERO-X]
        set(ff, 'position', [1250+sc2 570 560 420]);    % Time[PHASE] -> Time[ZERO-X]
        bode(-Gpd,opts); hold on; bode(Gpp,opts);  
        ylim([-360,0]); legend('Gpd', 'Gpp'); title('All Perturbations'); 
    end

    if(0)
        ff = figure(7); hold off;                  
        set(ff, 'position', [130+sc2 570 560 420]);    
        bode(-Gvd,opts); hold on; bode(-Gvd1,opts); bode(-Gvd2,opts); 
        bode(-Gvd3,opts); bode(-(Gvd1+Gvd2+Gvd3),opts); 
        legend('Gvd', 'Gvd1', 'Gvd2', 'Gvd3', 'Sum-1-2-3'); 
        ylim([-360,0]); title('Compare: Gvd');

        ff = figure(8); hold off;                   
        set(ff, 'position', [130+sc2 70 560 420]);   
        bode(Gvp,opts); hold on; bode(Gvp1,opts); bode(Gvp2,opts); 
        bode(Gvp3,opts); bode((Gvp1+Gvp2+Gvp3),opts); 
        legend('Gvp', 'Gvp1', 'Gvp2', 'Gvp3', 'Sum-1-2-3'); 
        ylim([-360,0]); title('Compare: Gvp');

        ff = figure(9); hold off;                 
        set(ff, 'position', [690+sc2 70 560 420]); 
        bode(-Gpd,opts); hold on; bode(-Gpd1,opts); bode(-Gpd2,opts); 
        bode(-Gpd3,opts); bode(-(Gpd1+Gpd2+Gpd3),opts); 
        legend('Gpd', 'Gpd1', 'Gpd2', 'Gpd3', 'Sum-1-2-3'); 
        ylim([-360,0]); title('Compare: Gpd');

        ff = figure(10); hold off;                   
        set(ff, 'position', [1250+sc2 70 560 420]); 
        bode(Gpp,opts); hold on; bode(Gpp1,opts); bode(Gpp2,opts);
        bode(Gpp3,opts); bode((Gpp1+Gpp2+Gpp3),opts); 
        legend('Gp', 'Gpp1', 'Gpp2', 'Gpp3', 'Sum-1-2-3'); 
        ylim([-360,0]); title('Compare: Gpp');
    end

    Gs = [Gvd, Gvp, Gpd, Gpp]; 
    GsName = {'Gvd', 'Gvp', 'Gpd', 'Gpp'}; 
    for i = 1:length(Gs)
    [Zc,Pc,~] = zpkdata(d2c(Gs(i),'tustin')); 
    Zc = Zc{1}; Pc = Pc{1};                             % unload from cell arrays           
    if(~isempty(Zc))
        inds = (real(Zc)>0);                            % All RHP zeros.
        Zfs = abs(Zc)/2/pi/1e3;                         % Frequencies. 
        if(any(inds) && ~isempty(rPair))                % If there's relevent zeros. 
            ind = find(min(Zfs(inds))==Zfs,1);          % Lowest freq. RHP zero.
            RHPzR = real(Zc(ind)); RHPzF = Zfs(ind);    % Save.
            fprintf([GsName{i} ', Real Part: %9.2f \n'],RHPzR);
            fprintf([GsName{i} ', Frequency: %9.3f [kHz]\n'], RHPzF); 
            disp('------------------------------------'); end; end
  
%     disp(['----------------------------- ' GsName{i} ': [S Domain]------------------------------']); 
%     disp('---------------[Zeros]--------------|-----------------[Poles]--------------'); 
%     disp('     [kHz]       [sig]        [jw]  |       [kHz]       [sig]        [jw]  '); 
%     for j = 1:length(real(Zc))  
%         fprintf('%10.3f  %10.0f  %10.0f  |  ', abs(Zc(j))/2/pi/1e3, real(Zc(j)), imag(Zc(j))); 
%         fprintf('%10.3f  %10.0f  %10.0f\n', abs(Pc(j))/2/pi/1e3, real(Pc(j)), imag(Pc(j))); end
%     disp('------------------------------------|--------------------------------------'); disp(' '); 
end
end
% -------------------------------------------------------------- 
        

%% ------------------ Freq. Sync. Loop Design ------------------- 
if(1)
    s = tf('s'); 
    z = tf('z', Pclk20); 

    Kpfd = (clks)/2/pi;                                 % Derived in notes. 
    a = 1000; DCOPclk = Pclk150;                        % a: # of clks where linearization is done. 
    Kdco = -1/a^2/DCOPclk;                              % DCOPclk: clk period of DCO in FPGA. 
    Gint = (Pclk20/2) * (1*z+1)/(z-1);                  % 1/s @ disc. sample time of Pclk20.
    [GintN, GintD] = tfdata(Gint, 'v'); 	

    K1 = 2^-28; K2 = 2^9; K3 = 2^-08; K6 = 2^-2;        % isoPhase #2: Dissertation ch7
    Gcpll = K1 * (K6/(1-z^-1) + K2/(1-(1-K3)*z^-1));    % PLL compensator.  
    [GcpllN, GcpllD] = tfdata(Gcpll, 'v'); 

%     K1 = 2^-19; K2 = 2^9; K3 = 2^-08; K6 = 2^-2;        % COMPEL 2020 DIGEST
%     Gcpll = K1 * (K6/(1-z^-1) + K2/(1-(1-K3)*z^-1));    % isoPhase: Dissertation ch7
%     [GcpllN, GcpllD] = tfdata(Gcpll, 'v'); 

    Gcplls = d2c(Gcpll, 'tustin');      [GcpllsN, GcpllsD] = tfdata(Gcplls); 
    Gints = d2c(Gint, 'tustin');        [GintsN, GintsD] = tfdata(Gints); 
    Gpps = d2c(Gpp, 'tustin');          [GppsN, GppsD] = tfdata(Gpps); 
    Gpds = d2c(Gpd, 'tustin');          [GpdsN, GpdsD] = tfdata(Gpds); 

    GfreqUC =(Gpps-1)*Kpfd*Kdco*2*pi*Gints;             % Uncompensated Loop, Tu. 
    Gfreq =GfreqUC*Gcplls;                              % Loop, T. 
    
    GfreqUCS = prescale(GfreqUC, {0.01,150e3});         % Avoid errors. 
    GfreqS = prescale(Gfreq, {0.01,150e3});    

    ff = figure(4); hold off;                           % Plot the system.
    set(ff, 'position', [0050 70 560 620]);
    bode(Gcplls,opts); hold on; 
    bode(Gpps-1,opts); bode(GfreqUCS,opts); margin(GfreqS); 
    ll = legend( 'Gcpll', 'Gpp-1', 'Freq. UnComp (Tu)', 'Freq. Loop (T)'); 
    ll.FontSize = 10;
    xlim([1e-4, 1.5e2]);
    
%     targ = -280;
%     [Ms, Ps, Ws] = bode(Gpps-1, logspace(1, 6, 10000));
%     Ps = Ps - ceil(max(Ps)/360)*360; 
%     ind = find(min(abs(Ps-targ))==abs(Ps-targ),1);
%     disp(['Plant Phase:   ' num2str(Ps(ind))]); 
%     disp(['Frequency:     ' num2str(Ws(ind)/2/pi)]); 
%     disp(['Phase Margin:  ' num2str(180+360+(-200)+Ps(ind))]); 
    
    % Publishable Graph
    if(1)
     	font = 'Times'; fntSz = 13; fStart = 100; fStop = 150e3; 
        fc = 4.03; phic = 62.4-180; yMin1 = -20; yMax1 = 50; yMin2 = -200; yMax2 = 0;
%      	font = 'Times'; fntSz = 13; fStart = 1; fStop = 150e3; 
%         fc = 0.0188; phic = 71.4-180; yMin1 = -20; yMax1 = 50; yMin2 = -200; yMax2 = 0;
        ff = figure(7); hold off; 
            set(ff, 'DefaultTextFontName', font, 'DefaultAxesFontName', font);
            set(ff, 'position', [0210 70 500 400]);
        freqs = logspace(0.5, 6, 10000); 
%         [Mag,Ang,Ws] = bode((Gpps-1),freqs,opts); 
        [Mag,Ang,Ws] = bode(GfreqUC,freqs,opts); 
            Ang = Ang - 360*0; 
            subplot(2,1,1); hold off; 
            semilogx(Ws/2/pi/1e3,squeeze(db(Mag)), 'lineWidth', 2);
            subplot(2,1,2); hold off; 
            semilogx(Ws/2/pi/1e3,squeeze(Ang), 'lineWidth', 2); 
        [Mag,Ang,Ws] = bode(Gcplls,freqs,opts); 
            Ang = Ang - 360*0; 
            subplot(2,1,1); hold on; 
            semilogx(Ws/2/pi/1e3,squeeze(db(Mag)), 'lineWidth', 2);
            subplot(2,1,2); hold on; 
            semilogx(Ws/2/pi/1e3,squeeze(Ang), 'lineWidth', 2); 
        [Mag,Ang,Ws] = bode(GfreqS,freqs,opts); 
            Ang = Ang - 360*0; 
            subplot(2,1,1); hold on; 
            hh = semilogx(Ws/2/pi/1e3,squeeze(db(Mag)), 'lineWidth', 2);
            subplot(2,1,2); hold on; 
            hh2 = semilogx(Ws/2/pi/1e3,squeeze(Ang), 'lineWidth', 2); 
        sp1 = subplot(2,1,1); 
            pos = get(sp1, 'Position');
            posnew = pos; posnew(2) = posnew(2)-0.06; set(sp1, 'Position', posnew);
            semilogx([min(Ws), max(Ws)]/2/pi/1e3,[0, 0],':k', 'lineWidth', 1); 
            semilogx(fc, 0, 'ok', 'lineWidth', 1); 
            semilogx([fc, fc],[0, yMin1],':k', 'lineWidth', 1); 
            xlim([fStart/1e3, fStop/1e3]); 
            ylim([yMin1, yMax1]); 
            ax = ancestor(hh, 'axes'); 
            set(gca, 'xtick', []); 
            set(gca, 'xticklabel', []); 
            yrule = ax.YAxis; yrule.FontSize = fntSz; 
            ylabel('Magnitude [dB]');
            text(1, -6.5, '4.03 kHz', 'fontSize', fntSz); subplot(2,1,2);
%             text(0.03, 6.0, '0.019 kHz', 'fontSize', fntSz); subplot(2,1,2); 
            semilogx([min(Ws), max(Ws)]/2/pi/1e3,[-180, -180],':k', 'lineWidth', 1); 
            semilogx([fc, fc],[yMax2, -180],':k', 'lineWidth', 1); 
            semilogx(fc, phic,'ok', 'lineWidth', 1); 
            xlim([fStart/1e3, fStop/1e3]); 
            ylim([yMin2, yMax2]); 
            ax = ancestor(hh2, 'axes'); 
            xrule = ax.XAxis; xrule.FontSize = fntSz; 
            yrule = ax.YAxis; yrule.FontSize = fntSz; 
            xlabel('Frequency [kHz]'); 
            ylabel('Phase [degree]'); 
%             ll = legend('{\itG_{zp}}-1', '{\itG_{cpll}}', '{\itT_{freq}}'); 
            ll = legend('{\it{T_{Ufreq}}}', '{\itG_{cpll}}', '{\itT_{freq}}'); 
            ll.FontSize = fntSz;
            text(4.7, -140, '62.4^\circ', 'fontSize', fntSz); 
%             text(0.025, -130, '71.4^\circ', 'fontSize', fntSz); 
    end
    
end
% -------------------------------------------------------------- 


%% ------------ Output Reg. Loop Design [Duty Only] ------------- 
if(1)
    s = tf('s'); 
    z = tf('z', Pclk20); 

    H2 = 2.5/5;                 	% Res. Divider: Vout/Vin = 2.5/5                      
    Kadc = 2^11/3.3;                % 11 bit ADC with 0-to-3.3V	 

    KmodD = 716/2048;   % Gain of single modulator block 
    if(b11Mod>1466)     KM3 = 0.3333;	KM2 = 0;    	KM1 = 0;
    elseif(b11Mod>999) 	KM3 = 0.0485;	KM2 = 0.3101;   KM1 = 0; 
    else                KM3 = 0.0284; 	KM2 = 0.0923;  	KM1 = 0.2333;  end

    Kd1 = -2^-14;   Kd2 = 2^6; 
    Kd3 = 2^-7;     Kd4 = 2^-5;     Kd5 = 0*2^-6; 
    Gcd = Kd1 * z^-1 * (1/(1-z^-1) + Kd2/(1-(1-Kd3)*z^-1)); 	% Gcd: duty to Vout compensator.   
    xPol = 1/(1-(1-Kd4-Kd5)*z^-1);                    
    Gcd = Gcd * xPol;       
    [GcdN, GcdD] = tfdata(Gcd, 'v'); 

    Kd1 = -2^-15;   Kd2 = 2^6; % COMPEL 2020 DIGEST             % isoVolt: Dissertation ch7
    Kd3 = 2^-7;     Kd4 = 2^-5;     Kd5 = 0*2^-6; 
    Gcd = Kd1 * z^-1 * (1/(1-z^-1) + Kd2/(1-(1-Kd3)*z^-1)); 	% Gcd: duty to Vout compensator.   
    xPol = 1/(1-(1-Kd4-Kd5)*z^-1);                    
    Gcd = Gcd * xPol;       
    [GcdN, GcdD] = tfdata(Gcd, 'v'); 

%     Kd1 = -2^-16;   Kd2 = 2^6;     
%     Kd3 = 2^-7;     Kd4 = 2^-5;     Kd5 = 0*2^-6; 
%     Gcd = Kd1 * z^-1 * (1/(1-z^-1) + Kd2/(1-(1-Kd3)*z^-1)); 	% Gcd: duty to Vout compensator.   
%     xPol = 1/(1-(1-Kd4-Kd5)*z^-1);                    
%     Gcd = Gcd * xPol;       
%     [GcdN, GcdD] = tfdata(Gcd, 'v'); 
    
%     Kd1 = -2^-19;   Kd2 = 2^10; 
%     Kd3 = 2^-5;     Kd4 = 2^-7;     Kd5 = 0*2^-6; 
%     Gcd = Kd1 * z^-1 * (1/(1-z^-1) + Kd2/(1-(1-Kd3)*z^-1)); 	% Gcd: duty to Vout compensator.   
%     xPol = 1/(1-(1-Kd4-Kd5)*z^-1);                    
%     Gcd = Gcd * xPol;       
%         [GcdN, GcdD] = tfdata(Gcd, 'v');    
%     %     figure(8); hold off; bode(Gcd, opts);
 	

    Gcd150 = d2d(Gcd, Ts, 'tustin');    [Gcd150N, Gcd150D] = tfdata(Gcd150); 
    Gcds = d2c(Gcd,'tustin');           [GcdsN, GcdsD] = tfdata(Gcds); 
    Gvds = d2c(Gvd, 'tustin');      
    Gvd1s = d2c(Gvd1, 'tustin');        [Gvd1sN, Gvd1sD] = tfdata(Gvd1s);  
    Gvd2s = d2c(Gvd2, 'tustin');        [Gvd2sN, Gvd2sD] = tfdata(Gvd2s); 
    Gvd3s = d2c(Gvd3, 'tustin');        [Gvd3sN, Gvd3sD] = tfdata(Gvd3s);  
    Gvps = d2c(Gvp, 'tustin');          [GvpsN, GvpsD] = tfdata(Gvps); 

    GvdAs = KM1*Gvd1s + KM2*Gvd2s + KM3*Gvd3s;          % All LVLs and Modulations. 

    GvLoopD = H2*Kadc*Gcds*Pclk150*1.0*GvdAs;               % Loop gain, T. 

    Tcl_vout = 1/Kadc/H2 * GvLoopD/(1+GvLoopD);
   	[mVout, mTime] = step(Tcl_vout,500*1e-6); 
    mTime = mTime - min(mTime); 
    mVout = ( mVout.*(-20) )/2; 
    
    ff = figure(5); hold off;                           % Plot the system.
    set(ff, 'position', [2610-2160 70 560 620]);
    bode(Gcds,opts); hold on;                           
    bode(GvdAs,opts); margin(GvLoopD); 
    ll = legend( 'Gcm', 'Gvm: 1+2+3',  'Loop: Tout'); 
    ll.FontSize = 10;
    xlim([1e-2, 1.5e2]);
    
    [pks, inds] = findpeaks(mVout); 
    pkTimes = mTime(inds); 
    pkPer = mean(diff(pkTimes(2:end))); 
    ringFreq = 1/pkPer; 
    
%     ff = figure(7); hold off;                           % Plot the system.
%     set(ff, 'position', [3260-3260 70 560 420]);
%     plot([0; 50+mTime*1e6], [mVout(1); mVout], 'lineWidth', 3); hold on; 
%     plot(pkTimes(2:end)*1e6+50, mVout(inds(2:end)), 'or', 'lineWidth', 2); 
%     title(sprintf('Ring Freq: %5.2f [kHz]', ringFreq/1e3)); 
%     ll.FontSize = 10;
    
%     targ = -280;
%     [Ms, Ps, Ws] = bode(Gvds, logspace(1, 6, 10000));
%     Ps = Ps - ceil(max(Ps)/360)*360; 
%     ind = find(min(abs(Ps-targ))==abs(Ps-targ),1);
%     disp(['Plant Phase:   ' num2str(Ps(ind))]); 
%     disp(['Frequency:     ' num2str(Ws(ind)/2/pi)]); 
%     disp(['Phase Margin:  ' num2str(180+360+(-200)+Ps(ind))]); 
    
    % Publishable graph. 
    if(1)
        font = 'Times'; fntSz = 13; fStart = 100; fStop = 150e3; 
        fc = 5.63; phic = 60.8-180; yMin1 = -40; yMax1 = 40; yMin2 = -360; yMax2 = 0; 
        ff = figure(7); hold off; 
            set(ff, 'DefaultTextFontName', font, 'DefaultAxesFontName', font);
            set(ff, 'position', [0210 70 500 400]);
        freqs = logspace(1, 6, 10000); 
%         [Mag,Ang,Ws] = bode(GvdAs*Pclk150/Ts,freqs,opts); 
        [Mag,Ang,Ws] = bode(H2*Kadc*Pclk150*GvdAs,freqs,opts); 
            Ang = Ang - 360*2; 
            subplot(2,1,1); hold off; 
                semilogx(Ws/2/pi/1e3,squeeze(db(Mag)), 'lineWidth', 2);
                subplot(2,1,2); hold off; 
                semilogx(Ws/2/pi/1e3,squeeze(Ang), 'lineWidth', 2); 
        [Mag,Ang,Ws] = bode(Gcds,freqs,opts); 
            Ang = Ang - 360*1; 
            subplot(2,1,1); hold on; 
                semilogx(Ws/2/pi/1e3,squeeze(db(Mag)), 'lineWidth', 2);
                subplot(2,1,2); hold on; 
                semilogx(Ws/2/pi/1e3,squeeze(Ang), 'lineWidth', 2); 
        [Mag,Ang,Ws] = bode(GvLoopD,freqs,opts); 
            Ang = Ang - 360*2; 
            subplot(2,1,1); hold on; 
                hh = semilogx(Ws/2/pi/1e3,squeeze(db(Mag)), 'lineWidth', 2);
                subplot(2,1,2); hold on; 
                hh2 = semilogx(Ws/2/pi/1e3,squeeze(Ang), 'lineWidth', 2); 
        sp1 = subplot(2,1,1); 
            pos = get(sp1, 'Position');
            posnew = pos; posnew(2) = posnew(2)-0.06; set(sp1, 'Position', posnew);
            semilogx([min(Ws), max(Ws)]/2/pi/1e3,[0, 0],':k', 'lineWidth', 1); 
            semilogx(fc, 0, 'ok', 'lineWidth', 1); 
            semilogx([fc, fc],[0, yMin1],':k', 'lineWidth', 1); 
            xlim([fStart/1e3, fStop/1e3]); 
            ylim([yMin1, yMax1]); 
            ax = ancestor(hh, 'axes'); 
            set(gca, 'xtick', []); 
            set(gca, 'xticklabel', []); 
            yrule = ax.YAxis; yrule.FontSize = fntSz; 
            ylabel('Magnitude [dB]');
            text(6.5, 7, '5.63 kHz', 'fontSize', fntSz); 
        subplot(2,1,2); 
            semilogx([min(Ws), max(Ws)]/2/pi/1e3,[-180, -180],':k', 'lineWidth', 1); 
            semilogx([fc, fc],[yMax2, -180],':k', 'lineWidth', 1); 
            semilogx(fc, phic,'ok', 'lineWidth', 1); 
            xlim([fStart/1e3, fStop/1e3]); 
            ylim([yMin2, yMax2]); 
            ax = ancestor(hh2, 'axes'); 
            xrule = ax.XAxis; xrule.FontSize = fntSz; 
            yrule = ax.YAxis; yrule.FontSize = fntSz; 
            xlabel('Frequency [kHz]'); 
            ylabel('Phase [degree]'); 
%             ll = legend('{\it\SigmaG_{vmn}}', '{\itG_{cm}}', '{\itT_{out}}'); 
            ll = legend('{\itT_{Uout}}', '{\itG_{cm}}', '{\itT_{out}}'); 
            ll.FontSize = fntSz;
            text(6.5, -80, '60.8^\circ', 'fontSize', fntSz); 
    end
    
%     ff = figure(6); hold off;                           % Plot the system.
%     set(ff, 'position', [3260 70 560 620]);
%     bode(1/Kadc/H2*GvLoopD/(1+GvLoopD)); 
%     ll = legend( 'T/1+T'); 
%     ll.FontSize = 10;
%     xlim([1e-2, 1.5e2]);
end
% -------------------------------------------------------------- 


% ------------------ Simulink Model Setup ----------------------
if(1)
    Az = top.As(:,:,1); 
    Ap1 = top.As(:,:,2); 
    Ap2 = top.As(:,:,3); 
    Ap3 = top.As(:,:,4); 
    An1 = top.As(:,:,5); 
    An2 = top.As(:,:,6); 
    An3 = top.As(:,:,7); 
    B = top.Bs(:,:,1); 
    initMod = b11Mod; 

    Gpd1s = d2c(Gpd1, 'tustin');    [Gpd1sN, Gpd1sD] = tfdata(Gpd1s);  
    Gpd2s = d2c(Gpd2, 'tustin');    [Gpd2sN, Gpd2sD] = tfdata(Gpd2s); 
    Gpd3s = d2c(Gpd3, 'tustin');    [Gpd3sN, Gpd3sD] = tfdata(Gpd3s);  

    % unload = load('dataSave/Gpp_3lvlModel_ph125_duty400.mat', 'Gpp');
    % Gpp_lvl3 = unload.Gpp; 
    % figure(4); hold on; bode(Gpp_lvl3,opts); 
    % [GppN_lvl3,GppD_lvl3] = tfdata(Gpp_lvl3,'v');
end
% -------------------------------------------------------------- 
% -------------------------------------------------------------- 
























%% ------------------ Running Fund. Phase of Step Response --------------------
if(0)
    % This code is useful for looking at the phase characteristics of a step
    % response... ie how does the real-time phase of the rectifier change with
    % a change in duty cycle. 
    if(0) 
        clear wVinv2 wVrec2 wx2 wIntM2 tIi time2 tI i
        pers = 30; sam = 1000;                  % periods to show. samples/period.
        time2 = linspace(0,Ts*pers,sam*pers);	% time vector to plot 
        prevIntM = 1; int = 1; x0 = Xss;        % initializations

        tdnom = round(td/Pclk150);                 % [25 350 125 25 350 125]
        tdnew = [0 400 100 0 400 100];          % new duty cycle operating point
        % tdnew = [13 350 137 13 350 137];        % new phase operating point

        tdnoms = repmat(tdnom, 1,30); 
        tdloc = [tdnom tdnom repmat(tdnew, 1,28)];
        % tdloc = [tdnom tdnom [0 375 125 0 400 100] repmat(tdnew, 1,27)];
        % tdloc = [tdnom tdnom [0 400 100 12 376 112] [0 400 100 0 400 100] tdnew tdnew tdnew tdnew tdnew tdnew];
        % tdloc = tdnoms; % Steady state. 

        teloc = cumsum(tdloc); teloc = teloc*Pclk150; tdloc = tdloc*Pclk150; 

        for i = 1:1:length(time2)
            int = find(time2(i)<teloc,1)-0;
            intM = mod(int,6);
                if(intM==0) intM = 6; end                   % mod(6,6) = 0; should be 6.
                if(isempty(int) && i==length(time2))        % final interval.
                        int = length(teloc); intM = 6; end      
                if(int==1) tI = time2(i);                   % when int-1 doesn't exist.
                else tI = time2(i) - teloc(int-1); end
            if(prevIntM ~= intM && i ~= 1) 
                x0 = expm(topA(:,:,prevIntM)*tdloc(prevInt))*x0+...
                    topA(:,:,prevIntM)\(expm(topA(:,:,prevIntM)*...
                    tdloc(prevInt))-eye(5))*topB(:,:,prevIntM)*U;
                end
            wx2(:,i) = expm(topA(:,:,intM)*tI)*x0+...
                topA(:,:,intM)\(expm(topA(:,:,intM)*tI)-eye(5))*topB(:,:,intM)*U;  
                if(intM==2) wVrec(i) = wx2(5,i); 
                elseif(intM==5) wVrec(i) = -wx2(5,i); 
                else wVrec(i) = 0; end
            wVinv2(i) = Us(intM); 
            prevIntM = intM; 
            prevInt = int; 
        end

        figure(100); 
        subplot(2,1,1);   
            hold off; 
            plot(time2*1e6, wx2(5,:), 'lineWidth', 2);
            ylabel('Vout'); 
        subplot(2,1,2);
            hold off; 
            plot(time2*1e6, wVinv2, 'lineWidth', 2); hold on;  
            plot(time2*1e6, wVrec, 'lineWidth', 2); 
        %     plot(time2*1e6, wx2(2,:), 'lineWidth', 2);
        %     plot(time2*1e6, wx2(4,:), 'lineWidth', 2);  
            plot(time2*1e6, wx2(1,:)*5, 'lineWidth', 2); hold on; 
            plot(time2*1e6, (wx2(1,:)-wx2(3,:))*10, 'lineWidth', 2);  
        %     legend('Inv.', 'Rec.', 'Cp', 'Cs', 'Ltx', 'Lrx'); 
            legend('Inv.', 'Rec.', 'Ltx*5', 'Lrx*10'); 
            xlabel('Time [{\mu}s]');
        % --------------------------------------------------------------


        % ----------------- Calculate Real Time Phases -----------------
        clear theta tper is Vrecl is_err Vr_err is_theta Vr_theta n4 n
        theta = flip((-180:0.1:179)*pi/180);                	% All theta values to be checked. 
        time3 = linspace(0,Ts,sam);                             % Basic timespan for one period.
        time4 = repmat(time3*f*2*pi, length(theta), 1);         % All time values for all theta instances. 
            time4 = time4 + theta';                             % Shift each row according to that theta. 
        wave = 15*sin(time4);                                   % All sin waves, each phase shifted. 
        nn = linspace(1,pers,200);                              % Check phase 4 times per period.
        shift = mod(nn,1) *2*pi;                                % Shift answers by this much.
        for n = 1:1:length(nn)                                  % Periods.
            tper = round((nn(n)-1)*sam)+1:1:round(nn(n)*sam);  	% Time indices for current period. 
            if(max(tper)>length(time2)) 
                tper = length(time2)-1000:1:length(time2); end
            is = (wx2(1,tper)-wx2(3,tper))*10;
            Vrecl = wVrec(tper); 
            Vinvl = wVinv2(tper); 
            Vcsl = wx2(4,tper); 
            Voutl(n) = max(wx2(5,tper)); 
            shift(n) = mod(nn(n),1) *pi; 
            for k = 1:1:length(theta)                           % Phases (theta).
                is_err(k) = sum(abs(is-wave(k,:)))/sam;         % Error relative to current theta. 
                Vr_err(k) = sum(abs(Vrecl-wave(k,:)))/sam;  
                Vcs_err(k) = sum(abs(Vcsl-wave(k,:)))/sam;   
                Vinv_err(k) = sum(abs(Vinvl-wave(k,:)))/sam;     
            end
            isThInd(n) = find(is_err == min(is_err),1); 	% Indices of min error (i.e.
            VrThInd(n) = find(Vr_err == min(Vr_err),1);  	% the phase of the fundamental).
            VcsThInd(n) = find(Vcs_err == min(Vcs_err),1);  	
            VinvThInd(n) = find(Vinv_err == min(Vinv_err),1);  
        end

        RecTh = theta(VrThInd)-theta(isThInd);                % The per-period fund. rect. input phase. 
        RecTh(RecTh>pi) = RecTh(RecTh>pi)-2*pi;         % Wrap Phase -180<ph<180.
        RecTh(RecTh<-pi) = RecTh(RecTh<-pi)+2*pi; 

        VcsTh = -theta(VcsThInd)+theta(VinvThInd)-pi/2+27.1*pi/180; 
    %     VcsTh = -theta(VcsThInd)+theta(VinvThInd); 
        VcsTh(VcsTh>pi) = VcsTh(VcsTh>pi)-2*pi; 
        VcsTh(VcsTh<-pi) = VcsTh(VcsTh<-pi)+2*pi; 

        % DEBUG: Check interval 'n' to ensure that the result looks right. 
        n = 17;     
            tper = round((nn(n)-1)*sam)+1:1:round(nn(n)*sam);   
            if(max(tper)>length(time2)) 
                tper = length(time2)-1000:1:length(time2); end
            is = (wx2(1,tper)-wx2(3,tper))*10;
            Vrecl = wVrec(tper); 
            Vcsl = wx2(4,tper)/3; 
            shift(n) = mod(nn(n),1) *pi;
        figure(101); hold off; 
        plot(time3, Vrecl, 'lineWidth', 2); hold on;  
        plot(time3, wave(VrThInd(n),:));
        plot(time3, is, 'lineWidth', 2);
        plot(time3, wave(isThInd(n),:)); 
        plot(time3, Vcsl, 'lineWidth', 2); 
        plot(time3, wave(VcsThInd(n),:)); 
        legend('Vr', 'Vr Fund.', 'is', 'is Fund.', 'Vcs', 'Vcs Fund.'); 
        disp('------------'); 
        disp([num2str(theta(VrThInd(n))*180/pi) ' | '...
            num2str(theta(isThInd(n))*180/pi) ' | '...
            num2str(RecTh(n)*180/pi) ' | '...
            num2str(theta(VcsThInd(n))*180/pi)]); 
        disp('Vr  |  is  |  Rec  |  Vcs');   

        figure(102); 
        subplot(2,1,1); hold off;  
            plot(Ts*nn*1e6, RecTh *180/pi, 'lineWidth', 2); hold on; 
            plot(Ts*nn*1e6, VcsTh*180/pi, 'lineWidth', 2);
            ylabel('Phase [deg.]');
            legend('Rectifier Phase', 'Vcs Phase', 'Vout Envelope'); 
        subplot(2,1,2); hold off; 
            plot(Ts*nn*1e6, Voutl, 'lineWidth', 2);
            legend('Vout Envelope'); 
    %         plot(time2*1e6, wVrec2, 'lineWidth', 2); hold on; 
    %         plot(time2*1e6, wx2(1,:)*5, 'lineWidth', 2);
    %         plot(time2*1e6, wx2(4,:), 'lineWidth', 2);
    %         legend('Vrec', 'Is', 'Vcs'); 
            xlabel('Time [{\mu}s]');
    end
end
% --------------------------------------------------------------












% --------------- Functions (Vrec Referenced) ------------------
% -------------------- 7 Level Intervals -----------------------
function [td,te,te2,topA,topB,inPair0,inPolD,inPolP,Vrs,Us] = ...
            modInts7(clk2,Pclk150,ph,M1,M2,M3,Az,Ap1,Ap2,Ap3,An1,An2,An3, ...
                                              Bz,Bp1,Bp2,Bp3,Bn1,Bn2,Bn3)

    FV = ones(size(Bz));
    FV(:,1) = -1;
                                          
    if(ph==M1 || ph==M2 || ph==M3)
        error(['Interval Definitions: edge times are invalid: ph=' num2str(ph)...
            ', M1=' num2str(M1) ', M2=' num2str(M2) ', M3=' num2str(M3)]); 
    elseif(abs(ph)<M1)
        if(ph>0)
            td = [ph,M1-ph,M2-M1,M3-M2,clk2-2*M3,M3-M2,M2-M1,M1]; td = [td td];
            topA = cat(3, Az,Az,Ap1,Ap2,Ap3,Ap2,Ap1,Az,  	Az, Az, An1, An2, An3, An2, An1, Az); 
            topB = cat(3,FV.*Bz,Bz,Bp1,Bp2,Bp3,Bp2,Bp1,Bz, 	...
                                        Bz,FV.*Bz,FV.*Bn1,FV.*Bn2,FV.*Bn3,FV.*Bn2,FV.*Bn1,FV.*Bz); 
            Vrs = [0,0,1,2,3,2,1,0,                         0,0,-1,-2,-3,-2,-1,0];
            Us = [-1,1,1,1,1,1,1,1,                         1,-1,-1,-1,-1,-1,-1,-1]; 
            inPair0 = [2,3,4,5,6,7; 3,4,5,6,7,8]; 
                inPair02 = inPair0 + 8; 
                inPair0 = [inPair0 inPair02]';
            inPolD = [1,1,1,0,0,0,                          1,1,1,0,0,0]'; 
            inPolP = [0,0,0,0,0,0,                          0,0,0,0,0,0]';
        else
            td = [ph,M1-ph,M2-M1,M3-M2,clk2-2*M3,M3-M2,M2-M1,M1]; td = [td td];
            topA = cat(3, Az,Az,Ap1,Ap2,Ap3,Ap2,Ap1,Az,     Az, Az, An1, An2, An3, An2, An1, Az); 
            topB = cat(3,FV.*Bz,Bz,Bp1,Bp2,Bp3,Bp2,Bp1,Bz, 	...
                                        Bz,FV.*Bz,FV.*Bn1,FV.*Bn2,FV.*Bn3,FV.*Bn2,FV.*Bn1,FV.*Bz); 
            Vrs = [0,0,1,2,3,2,1,0,                         0,0,-1,-2,-3,-2,-1,0];
            Us = [-1,1,1,1,1,1,1,1,                         1,-1,-1,-1,-1,-1,-1,-1]; 
            inPair0 = [2,3,4,5,6,7; 3,4,5,6,7,8]; 
                inPair02 = inPair0 + 8; 
                inPair0 = [inPair0 inPair02]';
            inPolD = [1,1,1,0,0,0,                          1,1,1,0,0,0]'; 
            inPolP = [0,0,0,0,0,0,                          0,0,0,0,0,0]';
        end
    elseif(abs(ph)<M2) 
        if(ph>0) 
            td = [M1,ph-M1,M2-ph,M3-M2,clk2-2*M3,M3-M2,M2-M1,M1]; td = [td td];
            topA = cat(3,Az,Ap1,Ap1,Ap2,Ap3,Ap2,Ap1,Az,     Az,An1,An1,An2,An3,An2,An1,Az); 
            topB = cat(3,FV.*Bz,FV.*Bp1,Bp1,Bp2,Bp3,Bp2,Bp1,Bz,   ...
                                    Bz,Bn1,FV.*Bn1,FV.*Bn2,FV.*Bn3,FV.*Bn2,FV.*Bn1,FV.*Bz); 
            Vrs = [0,1,1,2,3,2,1,0,                         0,-1,-1,-2,-3,-2,-1,0];
            Us = [-1,-1,1,1,1,1,1,1,                        1,1,-1,-1,-1,-1,-1,-1]; 
            inPair0 = [1,3,4,5,6,7; 2,4,5,6,7,8];             
                inPair02 = inPair0 + 8; 
                inPair0 = [inPair0 inPair02]';
            inPolD = [1,1,1,0,0,0,                          1,1,1,0,0,0]'; 
            inPolP = [0,0,0,0,0,0,                          0,0,0,0,0,0]'; 
        else
            td = [M1,ph-M1,M2-ph,M3-M2,clk2-2*M3,M3-M2,M2-M1,M1]; td = [td td];
            topA = cat(3,Az,Ap1,Ap1,Ap2,Ap3,Ap2,Ap1,Az,     Az,An1,An1,An2,An3,An2,An1,Az); 
            topB = cat(3,FV.*Bz,FV.*Bp1,Bp1,Bp2,Bp3,Bp2,Bp1,Bz,   ...
                                    Bz,Bn1,FV.*Bn1,FV.*Bn2,FV.*Bn3,FV.*Bn2,FV.*Bn1,FV.*Bz); 
            Vrs = [0,1,1,2,3,2,1,0,                         0,-1,-1,-2,-3,-2,-1,0];
            Us = [-1,-1,1,1,1,1,1,1,                        1,1,-1,-1,-1,-1,-1,-1]; 
            inPair0 = [1,3,4,5,6,7; 2,4,5,6,7,8];             
                inPair02 = inPair0 + 8; 
                inPair0 = [inPair0 inPair02]';
            inPolD = [1,1,1,0,0,0,                          1,1,1,0,0,0]'; 
            inPolP = [0,0,0,0,0,0,                          0,0,0,0,0,0]';
        end
    elseif(abs(ph)<M3)
        if(ph>0)
            td = [M1,M2-M1,ph-M2,M3-ph,clk2-2*M3,M3-M2,M2-M1,M1]; td = [td td];
            topA = cat(3,Az,Ap1,Ap2,Ap2,Ap3,Ap2,Ap1,Az,     Az,An1,An2,An2,An3,An2,An1,Az); 
            topB = cat(3,FV.*Bz,FV.*Bp1,FV.*Bp2,Bp2,Bp3,Bp2,Bp1,Bz, 	...
                                        Bz,Bn1,Bn2,FV.*Bn2,FV.*Bn3,FV.*Bn2,FV.*Bn1,FV.*Bz); 
            Vrs = [0,1,2,2,3,2,1,0,                         0,-1,-2,-2,-3,-2,-1,0];
            Us = [-1,-1,-1,1,1,1,1,1,                       1,1,1,-1,-1,-1,-1,-1]; 
            inPair0 = [1,2,4,5,6,7; 2,3,5,6,7,8];           
                inPair02 = inPair0 + 8; 
                inPair0 = [inPair0 inPair02]';
            inPolD = [1,1,1,0,0,0,                          1,1,1,0,0,0]'; 
            inPolP = [0,0,0,0,0,0,                          0,0,0,0,0,0]';
        else
            td = [M1,M2-M1,ph-M2,M3-ph,clk2-2*M3,M3-M2,M2-M1,M1]; td = [td td];
            topA = cat(3,Az,Ap1,Ap2,Ap2,Ap3,Ap2,Ap1,Az,     Az,An1,An2,An2,An3,An2,An1,Az); 
            topB = cat(3,FV.*Bz,FV.*Bp1,FV.*Bp2,Bp2,Bp3,Bp2,Bp1,Bz, 	...
                                        Bz,Bn1,Bn2,FV.*Bn2,FV.*Bn3,FV.*Bn2,FV.*Bn1,FV.*Bz); 
            Vrs = [0,1,2,2,3,2,1,0,                         0,-1,-2,-2,-3,-2,-1,0];
            Us = [-1,-1,-1,1,1,1,1,1,                       1,1,1,-1,-1,-1,-1,-1]; 
            inPair0 = [1,2,4,5,6,7; 2,3,5,6,7,8];           
                inPair02 = inPair0 + 8; 
                inPair0 = [inPair0 inPair02]';
            inPolD = [1,1,1,0,0,0,                          1,1,1,0,0,0]'; 
            inPolP = [0,0,0,0,0,0,                          0,0,0,0,0,0]';
        end
    elseif(abs(ph)>M3)
        if(ph>0)
            td = [M1,M2-M1,M3-M2,ph-M3,clk2-M3-ph,M3-M2,M2-M1,M1]; td = [td td];
            topA = cat(3,Az,Ap1,Ap2,Ap3,Ap3,Ap2,Ap1,Az,     Az,An1,An2,An3,An3,An2,An1,Az); 
            topB = cat(3,FV.*Bz,FV.*Bp1,FV.*Bp2,FV.*Bp3,Bp3,Bp2,Bp1,Bz,   ...
                                            Bz,Bn1,Bn2,Bn3,FV.*Bn3,FV.*Bn2,FV.*Bn1,FV.*Bz); 
            Vrs = [0,1,2,3,3,2,1,0,                         0,-1,-2,-3,-3,-2,-1,0];
            Us = [-1,-1,-1,-1,1,1,1,1,                      1,1,1,1,-1,-1,-1,-1]; 
            inPair0 = [1,2,3,5,6,7; 2,3,4,6,7,8];   
                inPair02 = inPair0 + 8; inPair0 = [inPair0 inPair02]';     
            inPolD = [1,1,1,0,0,0,                          1,1,1,0,0,0]'; 
            inPolP = [0,0,0,0,0,0,                          0,0,0,0,0,0]';
        else
            td = [M1,M2-M1,M3-M2,ph-M3,clk2-M3-ph,M3-M2,M2-M1,M1]; td = [td td];
            topA = cat(3,Az,Ap1,Ap2,Ap3,Ap3,Ap2,Ap1,Az,     Az,An1,An2,An3,An3,An2,An1,Az); 
            topB = cat(3,FV.*Bz,FV.*Bp1,FV.*Bp2,FV.*Bp3,Bp3,Bp2,Bp1,Bz,   ...
                                            Bz,Bn1,Bn2,Bn3,FV.*Bn3,FV.*Bn2,FV.*Bn1,FV.*Bz); 
            Vrs = [0,1,2,3,3,2,1,0,                         0,-1,-2,-3,-3,-2,-1,0];
            Us = [-1,-1,-1,-1,1,1,1,1,                      1,1,1,1,-1,-1,-1,-1]; 
            inPair0 = [1,2,3,5,6,7; 2,3,4,6,7,8];   
                inPair02 = inPair0 + 8; inPair0 = [inPair0 inPair02]';     
            inPolD = [1,1,1,0,0,0,                          1,1,1,0,0,0]'; 
            inPolP = [0,0,0,0,0,0,                          0,0,0,0,0,0]';
        end
    else
        error('Interval Definitions: unexpected if/else evaluation!!'); 
    end
    td = td*Pclk150; te = cumsum(td); te2 = [0 te]; 
end
% --------------------------------------------------------------

% --------------------------------------------------------------
function [M1, M2, M3] = fundToMs(fund)

    M1_list = [250,249,249,249,248,248,248,247,247,247,246,246,246,246,245,245,245,244,244,244,243,243,243,242,242,242,241,241,241,240,240,240,239,239,239,239,238,238,238,237,237,237,236,236,236,235,235,235,234,234,234,233,233,233,232,232,232,232,231,231,231,230,230,230,229,229,229,228,228,228,227,227,227,226,226,226,225,225,225,224,224,224,224,223,223,223,222,222,222,221,221,221,220,220,220,219,219,219,218,218,218,217,217,217,217,216,216,216,215,215,215,214,214,214,213,213,213,212,212,212,211,211,211,210,210,210,209,209,209,209,208,208,208,207,207,207,206,206,206,205,205,205,204,204,204,203,203,203,202,202,202,202,201,201,201,200,200,200,199,199,199,198,198,198,197,197,197,196,196,196,195,195,195,195,194,194,194,193,193,193,192,192,192,191,191,191,190,190,190,189,189,189,188,188,188,187,187,187,187,186,186,186,185,185,185,184,184,184,183,183,183,182,182,182,181,181,181,180,180,180,180,179,179,179,178,178,178,177,177,177,176,176,176,175,175,175,174,174,174,173,173,173,173,172,172,172,171,171,171,170,170,170,169,169,169,168,168,168,167,167,167,166,166,166,165,165,165,165,164,164,164,163,163,163,162,162,162,161,161,161,160,160,160,159,159,159,158,158,158,158,157,157,157,156,156,156,155,155,155,154,154,154,153,153,153,152,152,152,151,151,151,151,150,150,150,149,149,149,148,148,148,147,147,147,146,146,146,145,145,145,144,144,144,143,143,143,143,142,142,142,141,141,141,140,140,140,139,139,139,138,138,138,137,137,137,136,136,136,136,135,135,135,134,134,134,133,133,133,132,132,132,131,131,131,130,130,130,129,129,129,128,128,128,128,127,127,127,126,126,126,125,125,125,124,124,124,123,123,123,122,122,122,122,121,121,121,120,120,119,119,119,118,118,118,117,117,116,116,116,115,115,115,114,114,114,113,113,113,112,112,112,111,111,110,110,110,109,109,109,108,108,107,107,107,106,106,106,105,105,105,104,104,104,103,103,102,102,102,101,101,101,100,100,100,99,99,99,98,98,97,97,97,96,96,96,95,95,95,94,94,94,93,93,92,92,92,91,91,91,90,90,90,89,89,89,88,88,87,87,87,86,86,86,85,85,85,84,84,84,83,83,82,82,82,81,81,81,80,80,80,79,79,78,78,78,77,77,77,76,76,76,75,75,75,74,74,73,73,73,72,72,72,71,71,71,70,70,70,69,69,68,68,68,67,67,67,66,66,66,65,65,65,64,64,64,63,63,62,62,62,61,61,61,60,60,60,59,59,58,58,58,57,57,57,56,56,56,55,55,55,55,55,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,48,48,48,48,48,48,48,48,48,48,48,48,48,48,48,48,48,48,48,48,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,44,44,44,44,44,44,44,44,44,44,44,44,44,44,44,44,44,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,41,41,41,41,41,41,41,41,41,41,41,41,41,41,41,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,38,38,38,38,38,38,38,38,38,38,38,38,38,38,38,37,37,37,37,37,37,37,37,37,37,37,37,37,37,36,36,36,36,36,36,36,36,36,36,36,36,36,35,35,35,35,35,35,35,35,35,35,35,35,35,35,34,34,34,34,34,34,34,34,34,34,34,34,34,34,33,33,33,33,33,33,33,33,33,33,33,33,33,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,4,4,4,4,4,4,4,4,4,4,4,4,4,4];
    M2_list = [250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,249,249,249,248,248,248,247,247,247,247,246,246,246,245,245,245,244,244,244,243,243,243,242,242,242,241,241,241,241,240,240,240,239,239,239,238,238,238,237,237,237,236,236,236,236,235,235,235,234,234,234,233,233,233,232,232,232,231,231,231,230,230,230,230,229,229,229,228,228,228,227,227,227,226,226,226,225,225,225,224,224,224,224,223,223,223,222,222,222,221,221,221,220,220,220,219,219,219,218,218,218,218,217,217,217,216,216,216,215,215,215,214,214,214,213,213,213,213,212,212,212,211,211,211,210,210,210,209,209,209,208,208,208,207,207,207,207,206,206,206,205,205,205,204,204,204,203,203,203,202,202,202,201,201,201,201,200,200,200,199,199,199,198,198,198,197,197,197,196,196,196,195,195,195,195,194,194,194,193,193,193,192,192,192,191,191,191,190,190,190,190,189,189,189,188,188,188,187,187,187,187,187,186,186,186,186,185,185,184,184,184,184,183,183,183,182,182,182,182,181,181,180,180,180,180,179,179,178,178,178,178,177,177,177,176,176,176,176,175,175,175,174,174,174,173,173,173,172,172,172,171,171,171,171,170,170,170,169,169,169,168,168,168,167,167,167,167,166,166,166,165,165,165,164,164,164,164,163,163,163,162,162,162,161,161,161,161,160,160,160,159,159,159,159,158,158,157,157,157,156,156,156,156,155,155,155,155,154,154,154,153,153,153,153,152,152,151,151,151,151,150,150,150,149,149,149,149,148,148,148,147,147,147,147,146,146,146,145,145,145,144,144,144,144,143,143,143,142,142,142,141,141,141,141,140,140,140,139,139,139,139,138,138,138,137,137,137,136,136,136,136,135,135,135,134,134,134,134,133,133,133,132,132,132,132,131,131,131,130,130,130,129,129,129,129,128,128,128,127,127,127,127,126,126,126,125,125,125,125,124,124,124,124,123,123,123,122,122,122,122,121,121,121,120,120,120,119,119,119,119,118,118,118,117,117,117,117,116,116,116,115,115,115,115,114,114,114,113,113,113,113,112,112,112,111,111,111,111,110,110,110,110,109,109,109,108,108,108,108,107,107,107,106,106,106,106,105,105,105,104,104,104,104,103,103,103,103,103,103,103,103,103,103,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,101,101,101,101,101,101,101,101,101,101,101,101,101,101,101,101,101,101,101,101,101,101,101,101,101,101,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,98,98,98,98,98,98,98,98,98,98,98,98,98,98,98,98,98,97,97,97,97,97,97,97,97,97,97,97,97,97,97,97,97,97,96,96,96,96,96,96,96,96,96,96,96,96,96,96,96,95,95,95,95,95,95,95,95,95,95,95,95,95,95,95,94,94,94,94,94,94,94,94,94,94,94,94,94,94,93,93,93,93,93,93,93,93,93,93,93,93,93,92,92,92,92,92,92,92,92,92,92,92,92,92,91,91,91,91,91,91,91,91,91,91,91,91,91,90,90,90,90,90,90,90,90,90,90,90,90,89,89,89,89,89,89,89,89,89,89,89,89,88,88,88,88,88,88,88,88,88,88,88,88,87,87,87,87,87,87,87,87,87,87,87,87,86,86,86,86,86,86,86,86,86,86,86,85,85,85,85,85,85,85,85,85,85,85,84,84,84,84,84,84,84,84,84,84,84,83,83,83,83,83,83,83,83,83,83,83,82,82,82,82,82,82,82,82,82,82,82,81,81,81,81,81,81,81,81,81,81,81,80,80,80,80,80,80,80,80,80,80,79,79,79,79,79,79,79,79,79,79,78,78,78,78,78,78,78,78,78,78,78,77,77,77,77,77,77,77,77,77,77,76,76,76,76,76,76,76,76,76,76,75,75,75,75,75,75,75,75,75,75,74,74,74,74,74,74,74,74,74,74,73,73,73,73,73,73,73,73,73,73,72,72,72,72,72,72,72,72,72,72,71,71,71,71,71,71,71,71,71,71,70,70,70,70,70,70,70,70,70,69,69,69,69,69,69,69,69,69,69,68,68,68,68,68,68,68,68,68,67,67,67,67,67,67,67,67,67,67,66,66,66,66,66,66,66,66,66,65,65,65,65,65,65,65,65,65,64,64,64,64,64,64,64,64,64,64,63,63,63,63,63,63,63,63,63,62,62,62,62,62,62,62,62,62,61,61,61,61,61,61,61,61,61,61,60,60,60,60,60,60,60,60,60,59,59,59,59,59,59,59,59,59,58,58,58,58,58,58,58,58,58,57,57,57,57,57,57,57,57,57,56,56,56,56,56,56,56,56,56,55,55,55,55,55,55,55,55,55,54,54,54,54,54,54,54,54,54,53,53,53,53,53,53,53,53,53,52,52,52,52,52,52,52,52,52,51,51,51,51,51,51,51,51,51,50,50,50,50,50,50,50,50,50,49,49,49,49,49,49,49,49,49,48,48,48,48,48,48,48,48,48,47,47,47,47,47,47,47,47,46,46,46,46,46,46,46,46,46,45,45,45,45,45,45,45,45,45,44,44,44,44,44,44,44,44,44,43,43,43,43,43,43,43,43,43,42,42,42,42,42,42,42,42,41,41,41,41,41,41,41,41,41,40,40,40,40,40,40,40,40,40,39,39,39,39,39,39,39,39,38,38,38,38,38,38,38,38,38,37,37,37,37,37,37,37,37,37,36,36,36,36,36,36,36,36,35,35,35,35,35,35,35,35,35,34,34,34,34,34,34,34,34,34,33,33,33,33,33,33,33,33,32,32,32,32,32,32,32,32,32,31,31,31,31,31,31,31,31,30,30,30,30,30,30,30,30,30,29,29,29,29,29,29,29,29,29,28,28,28,28,28,28,28,28,27,27,27,27,27,27,27,27,27,26,26,26,26,26,26,26,26,25,25,25,25,25,25,25,25,25,24,24,24,24,24,24,24,24,23,23,23,23,23,23,23,23,23,22,22,22,22,22,22,22,22,21,21,21,21,21,21,21,21,21,20,20,20,20,20,20,20,20,19,19,19,19,19,19,19,19,19,18,18,18,18,18,18,18,18,17,17,17,17,17,17,17,17,17,16,16,16,16,16,16,16,16,15,15,15,15,15,15,15,15,14,14,14,14,14,14,14,14,14,13,13,13,13,13,13,13,13,12,12,12,12,12,12];
    M3_list = [250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,249,249,249,249,248,248,248,248,247,247,247,246,246,246,246,245,245,245,245,244,244,244,243,243,243,243,242,242,242,241,241,241,241,240,240,240,240,239,239,239,238,238,238,238,237,237,237,237,236,236,236,235,235,235,235,234,234,234,234,233,233,233,232,232,232,232,231,231,231,230,230,230,230,229,229,229,229,228,228,228,227,227,227,227,226,226,226,226,225,225,225,224,224,224,224,223,223,223,223,222,222,222,221,221,221,221,220,220,220,219,219,219,219,218,218,218,218,217,217,217,216,216,216,216,215,215,215,215,214,214,214,213,213,213,213,212,212,212,212,211,211,211,210,210,210,210,209,209,209,208,208,208,208,207,207,207,207,206,206,206,205,205,205,205,204,204,204,204,203,203,203,202,202,202,202,201,201,201,201,200,200,200,199,199,199,199,198,198,198,197,197,197,197,196,196,196,196,195,195,195,195,195,195,194,194,194,194,193,193,193,192,192,192,192,191,191,191,191,190,190,190,189,189,189,189,188,188,188,187,187,187,187,186,186,186,186,185,185,185,185,184,184,184,184,183,183,183,182,182,182,182,182,181,181,181,180,180,180,180,179,179,179,179,178,178,178,178,177,177,177,177,176,176,176,176,175,175,175,175,174,174,174,174,173,173,173,173,172,172,172,172,171,171,171,171,170,170,170,170,169,169,169,169,168,168,168,168,167,167,167,167,166,166,166,166,165,165,165,165,165,164,164,164,163,163,163,163,163,162,162,162,162,161,161,161,161,160,160,160,160,159,159,159,159,158,158,158,158,157,157,157,157,156,156,156,156,156,155,155,155,155,154,154,154,154,153,153,153,153,152,152,152,152,152,151,151,151,151,150,150,150,150,149,149,149,149,148,148,148,148,148,147,147,147,147,146,146,146,146,145,145,145,145,144,144,144,144,144,143,143,143,143,142,142,142,142,142,141,141,141,141,140,140,140,140,139,139,139,139,139,138,138,138,138,137,137,137,137,136,136,136,136,136,135,135,135,135,134,134,134,134,134,133,133,133,133,132,132,132,132,132,131,131,131,131,130,130,130,130,129,129,129,129,129,128,128,128,128,128,127,127,127,127,126,126,126,126,126,125,125,125,125,124,124,124,124,123,123,123,123,123,122,122,122,122,122,121,121,121,121,120,120,120,120,120,119,119,119,119,118,118,118,118,118,117,117,117,117,117,116,116,116,116,115,115,115,115,115,114,114,114,114,113,113,113,113,113,112,112,112,112,112,111,111,111,111,110,110,110,110,110,109,109,109,109,108,108,108,108,108,107,107,107,107,107,106,106,106,106,105,105,105,105,105,104,104,104,104,104,103,103,103,103,103,102,102,102,102,101,101,101,101,101,100,100,100,100,100,99,99,99,99,98,98,98,98,98,97,97,97,97,97,96,96,96,96,96,95,95,95,95,94,94,94,94,94,93,93,93,93,93,92,92,92,92,91,91,91,91,91,90,90,90,90,90,89,89,89,89,89,88,88,88,88,88,87,87,87,87,87,86,86,86,86,85,85,85,85,85,84,84,84,84,84,83,83,83,83,83,82,82,82,82,82,81,81,81,81,81,80,80,80,80,79,79,79,79,79,78,78,78,78,78,77,77,77,77,77,76,76,76,76,76,75,75,75,75,75,74,74,74,74,74,73,73,73,73,72,72,72,72,72,71,71,71,71,71,70,70,70,70,70,69,69,69,69,69,68,68,68,68,68,67,67,67,67,67,66,66,66,66,66,65,65,65,65,65,64,64,64,64,64,63,63,63,63,63,62,62,62,62,61,61,61,61,61,60,60,60,60,60,59,59,59,59,59,58,58,58,58,58,57,57,57,57,57,56,56,56,56,56,55,55,55,55,55,54,54,54,54,54,53,53,53,53,53,52,52,52,52,52,51,51,51,51,51,50,50,50,50,50,49,49,49,49,49,48,48,48,48,48,47,47,47,47,47,46,46,46,46,46,45,45,45,45,45,44,44,44,44,44,43,43,43,43,43,42,42,42,42,42,41,41,41,41,41,40,40,40,40,40,39,39,39,39,39,38,38,38,38,38,37,37,37,37,37,36,36,36,36,36,35,35,35,35,35,34,34,34,34,34,33,33,33,33,33,32,32,32,32,32,31,31,31,31,31,30,30,30,30,30,29,29,29,29,29,28,28,28,28,28,27,27,27,27,27,26,26,26,26,26,25,25,25,25,25,24,24,24,24,24,23,23,23,23,23,22,22,22,22,22,21,21,21,21,21,21,20,20,20];
    
    nA_list = [2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.51, 2.51, 2.52, 2.52, 2.53, 2.53, 2.54, 2.54, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.56, 2.56, 2.57, 2.57, 2.58, 2.58, 2.59, 2.59, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.61, 2.61, 2.62, 2.62, 2.63, 2.63, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.65, 2.65, 2.66, 2.66, 2.67, 2.67, 2.68, 2.68, 2.69, 2.69, 2.69, 2.69, 2.69, 2.69, 2.69, 2.69, 2.69, 2.69, 2.69, 2.70, 2.70, 2.70, 2.71, 2.71, 2.72, 2.72, 2.73, 2.73, 2.74, 2.74, 2.74, 2.74, 2.74, 2.74, 2.74, 2.74, 2.74, 2.75, 2.75, 2.76, 2.76, 2.77, 2.77, 2.78, 2.78, 2.79, 2.79, 2.79, 2.79, 2.79, 2.79, 2.79, 2.80, 2.80, 2.80, 2.81, 2.81, 2.82, 2.82, 2.83, 2.83, 2.83, 2.83, 2.83, 2.84, 2.84, 2.85, 2.85, 2.86, 2.86, 2.87, 2.87, 2.88, 2.88, 2.88, 2.88, 2.88, 2.89, 2.89, 2.90, 2.90, 2.90, 2.91, 2.91, 2.92, 2.92, 2.93, 2.93, 2.93, 2.93, 2.94, 2.94, 2.95, 2.95, 2.96, 2.96, 2.97, 2.97, 2.98, 2.98, 2.98, 2.99, 2.99, 3.00, 3.00, 3.00, 3.01, 3.01, 3.02, 3.02, 3.03, 3.03, 3.04, 3.04, 3.05, 3.05, 3.06, 3.06, 3.07, 3.08, 3.08, 3.09, 3.09, 3.10, 3.10, 3.10, 3.11, 3.11, 3.12, 3.13, 3.13, 3.14, 3.14, 3.15, 3.15, 3.16, 3.16, 3.17, 3.18, 3.18, 3.19, 3.19, 3.20, 3.20, 3.20, 3.21, 3.22, 3.23, 3.23, 3.24, 3.24, 3.25, 3.26, 3.27, 3.28, 3.28, 3.29, 3.29, 3.30, 3.30, 3.31, 3.32, 3.32, 3.33, 3.33, 3.34, 3.34, 3.35, 3.36, 3.37, 3.38, 3.38, 3.39, 3.40, 3.40, 3.41, 3.42, 3.43, 3.43, 3.44, 3.45, 3.46, 3.47, 3.47, 3.48, 3.48, 3.49, 3.50, 3.51, 3.52, 3.52, 3.53, 3.54, 3.55, 3.56, 3.57, 3.57, 3.58, 3.59, 3.59, 3.60, 3.61, 3.62, 3.63, 3.64, 3.65, 3.66, 3.67, 3.68, 3.69, 3.69, 3.70, 3.71, 3.72, 3.73, 3.74, 3.75, 3.76, 3.77, 3.78, 3.79, 3.80, 3.81, 3.82, 3.83, 3.84, 3.85, 3.86, 3.87, 3.88, 3.89, 3.89, 3.91, 3.92, 3.93, 3.94, 3.95, 3.97, 3.98, 3.99, 4.00, 4.01, 4.02, 4.03, 4.05, 4.06, 4.07, 4.08, 4.09, 4.11, 4.12, 4.13, 4.14, 4.16, 4.17, 4.18, 4.19, 4.20, 4.21, 4.23, 4.24, 4.26, 4.27, 4.28, 4.29, 4.31, 4.32, 4.34, 4.35, 4.36, 4.38, 4.39, 4.40, 4.42, 4.43, 4.45, 4.46, 4.48, 4.49, 4.51, 4.52, 4.54, 4.56, 4.57, 4.59, 4.60, 4.62, 4.63, 4.65, 4.67, 4.68, 4.70, 4.71, 4.73, 4.75, 4.77, 4.79, 4.80, 4.82, 4.84, 4.86, 4.88, 4.89, 4.91, 4.93, 4.95, 4.97, 4.99, 5.01, 5.03, 5.05, 5.07, 5.09, 5.11, 5.13, 5.15, 5.17, 5.19, 5.21, 5.24, 5.26, 5.28, 5.30, 5.32, 5.35, 5.37, 5.39, 5.41, 5.43, 5.46, 5.48, 5.51, 5.53, 5.56, 5.58, 5.61, 5.64, 5.67, 5.69, 5.71, 5.74, 5.77, 5.79, 5.82, 5.85, 5.88, 5.90, 5.93, 5.97, 5.99, 6.02, 6.05, 6.08, 6.11, 6.14, 6.17, 6.20, 6.24, 6.27, 6.30, 6.33, 6.37, 6.40, 6.43, 6.47, 6.51, 6.54, 6.58, 6.61, 6.65, 6.69, 6.73, 6.77, 6.80, 6.84, 6.88, 6.92, 6.97, 7.00, 7.04, 7.08, 7.12, 7.17, 7.21, 7.25, 7.29, 7.32, 7.35, 7.38, 7.40, 7.43, 7.45, 7.46, 7.47, 7.47, 7.48, 7.48, 7.49, 7.49, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.51, 7.51, 7.52, 7.52, 7.53, 7.53, 7.54, 7.54, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.56, 7.56, 7.57, 7.57, 7.57, 7.58, 7.58, 7.59, 7.59, 7.59, 7.59, 7.59, 7.59, 7.59, 7.59, 7.59, 7.60, 7.60, 7.61, 7.61, 7.62, 7.62, 7.63, 7.63, 7.64, 7.64, 7.64, 7.64, 7.64, 7.64, 7.64, 7.65, 7.65, 7.66, 7.66, 7.67, 7.67, 7.67, 7.68, 7.68, 7.69, 7.69, 7.69, 7.69, 7.70, 7.70, 7.71, 7.71, 7.72, 7.72, 7.73, 7.73, 7.74, 7.74, 7.74, 7.74, 7.75, 7.75, 7.76, 7.76, 7.77, 7.77, 7.77, 7.78, 7.78, 7.79, 7.79, 7.80, 7.80, 7.81, 7.81, 7.82, 7.82, 7.83, 7.83, 7.84, 7.84, 7.85, 7.85, 7.86, 7.86, 7.87, 7.87, 7.87, 7.88, 7.89, 7.89, 7.90, 7.90, 7.91, 7.91, 7.92, 7.92, 7.93, 7.94, 7.94, 7.95, 7.95, 7.96, 7.96, 7.97, 7.97, 7.98, 7.99, 7.99, 8.00, 8.00, 8.01, 8.01, 8.02, 8.03, 8.04, 8.04, 8.05, 8.05, 8.06, 8.07, 8.07, 8.08, 8.09, 8.09, 8.10, 8.10, 8.11, 8.12, 8.13, 8.14, 8.14, 8.15, 8.15, 8.16, 8.17, 8.18, 8.18, 8.19, 8.19, 8.20, 8.21, 8.22, 8.23, 8.23, 8.24, 8.24, 8.25, 8.26, 8.27, 8.28, 8.28, 8.29, 8.30, 8.31, 8.32, 8.33, 8.33, 8.34, 8.35, 8.36, 8.37, 8.37, 8.38, 8.38, 8.39, 8.40, 8.41, 8.42, 8.43, 8.44, 8.45, 8.46, 8.47, 8.47, 8.48, 8.49, 8.50, 8.51, 8.52, 8.52, 8.53, 8.54, 8.55, 8.56, 8.57, 8.58, 8.59, 8.60, 8.61, 8.62, 8.63, 8.64, 8.65, 8.66, 8.67, 8.67, 8.68, 8.69, 8.70, 8.71, 8.72, 8.73, 8.74, 8.75, 8.76, 8.77, 8.78, 8.79, 8.80, 8.81, 8.82, 8.83, 8.84, 8.85, 8.87, 8.87, 8.88, 8.89, 8.91, 8.92, 8.93, 8.94, 8.95, 8.97, 8.97, 8.98, 9.00, 9.01, 9.02, 9.03, 9.04, 9.06, 9.07, 9.07, 9.08, 9.10, 9.11, 9.12, 9.13, 9.14, 9.16, 9.17, 9.18, 9.19, 9.21, 9.22, 9.23, 9.25, 9.26, 9.27, 9.28, 9.30, 9.31, 9.32, 9.33, 9.35, 9.36, 9.37, 9.38, 9.40, 9.41, 9.42, 9.44, 9.45, 9.46, 9.48, 9.49, 9.50, 9.52, 9.53, 9.55, 9.56, 9.57, 9.59, 9.60, 9.62, 9.63, 9.65, 9.66, 9.67, 9.69, 9.70, 9.72, 9.73, 9.75, 9.76, 9.77, 9.79, 9.81, 9.82, 9.84, 9.86, 9.87, 9.89, 9.90, 9.92, 9.94, 9.95, 9.96, 9.98,10.00,10.01,10.03,10.05,10.06,10.08,10.09,10.11,10.13,10.14,10.16,10.18,10.19,10.21,10.23,10.25,10.26,10.28,10.30,10.32,10.34,10.36,10.37,10.39,10.41,10.43,10.45,10.46,10.48,10.50,10.52,10.54,10.56,10.57,10.59,10.61,10.63,10.65,10.67,10.69,10.71,10.73,10.75,10.77,10.79,10.81,10.83,10.85,10.87,10.89,10.91,10.93,10.96,10.97,10.99,11.02,11.04,11.06,11.08,11.10,11.12,11.15,11.16,11.19,11.21,11.24,11.26,11.28,11.31,11.33,11.35,11.38,11.40,11.43,11.45,11.47,11.50,11.52,11.55,11.57,11.59,11.62,11.64,11.66,11.69,11.71,11.74,11.76,11.78,11.81,11.83,11.85,11.88,11.91,11.94,11.96,11.99,12.01,12.04,12.07,12.10,12.13,12.15,12.18,12.21,12.24,12.26,12.29,12.32,12.35,12.37,12.39,12.41,12.43,12.44,12.45,12.46,12.47,12.48,12.48,12.49,12.49,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.51,12.51,12.52,12.52,12.53,12.53,12.54,12.54,12.55,12.55,12.55,12.55,12.55,12.55,12.55,12.55,12.55,12.55,12.55,12.55,12.55,12.56,12.56,12.57,12.57,12.58,12.58,12.59,12.59,12.59,12.59,12.59,12.59,12.59,12.60,12.60,12.61,12.61,12.62,12.62,12.63,12.63,12.64,12.64,12.64,12.64,12.65,12.65,12.65,12.66,12.66,12.67,12.67,12.68,12.68,12.69,12.69,12.69,12.70,12.70,12.71,12.71,12.72,12.72,12.73,12.73,12.74,12.74,12.75,12.75,12.75,12.76,12.76,12.77,12.77,12.78,12.79,12.79,12.80,12.80,12.81,12.81,12.82,12.82,12.83,12.84,12.85,12.85,12.85,12.86,12.86,12.87,12.88,12.89,12.89,12.90,12.90,12.91,12.91,12.92,12.93,12.94,12.94,12.95,12.95,12.95,12.96,12.97,12.98,12.99,12.99,13.00,13.00,13.01,13.02,13.03,13.04,13.04,13.05,13.05,13.05,13.06,13.07,13.08,13.09,13.09,13.10,13.10,13.11,13.12,13.13,13.14,13.14,13.15,13.15,13.16,13.17,13.18,13.19,13.19,13.20,13.21,13.22,13.23,13.24,13.24,13.25,13.26,13.27,13.28,13.28,13.29,13.30,13.31,13.32,13.33,13.33,13.34,13.35,13.36,13.37,13.38,13.38,13.39,13.40,13.41,13.42,13.43,13.44,13.44,13.45,13.46,13.47,13.48,13.49,13.50,13.51,13.52,13.53,13.54,13.55,13.56,13.57,13.58,13.59,13.60,13.61,13.62,13.63,13.64,13.64,13.65,13.66,13.67,13.68,13.69,13.70,13.72,13.73,13.74,13.74,13.75,13.77,13.78,13.79,13.80,13.81,13.82,13.83,13.84,13.84,13.86,13.87,13.88,13.89,13.90,13.92,13.93,13.94,13.94,13.96,13.97,13.98,13.99,14.00,14.01,14.02,14.03,14.04,14.06,14.07,14.08,14.09,14.11,14.12,14.13,14.14,14.15,14.16,14.17,14.19,14.20,14.21,14.22,14.24,14.25,14.26,14.27,14.28,14.30,14.31,14.32,14.34,14.35,14.36,14.37,14.39,14.40,14.41,14.43,14.44,14.45,14.46,14.47,14.49,14.50,14.51,14.53,14.54,14.55,14.56,14.58,14.59,14.61,14.62,14.64,14.64,14.66,14.67,14.69,14.70,14.72,14.73,14.74,14.76,14.77,14.79,14.80,14.82,14.83,14.84,14.86,14.87,14.89,14.90,14.92,14.94,14.95,14.96,14.98,14.99,15.01,15.02,15.04,15.05,15.06,15.08,15.09,15.11,15.13,15.14,15.15,15.17,15.19,15.20,15.22,15.23,15.25,15.27,15.28,15.30,15.32,15.33,15.34,15.36,15.38,15.39,15.41,15.43,15.44,15.46,15.48,15.50,15.52,15.53,15.55,15.57,15.58,15.60,15.62,15.63,15.65,15.66,15.68,15.70,15.72,15.73,15.75,15.77,15.79,15.81,15.83,15.84,15.86,15.88,15.90,15.91,15.93,15.95,15.97,15.99,16.00,16.02,16.04,16.06,16.08,16.10,16.12,16.13,16.15,16.17,16.19,16.21,16.23,16.25,16.27,16.29,16.31,16.33,16.35,16.37,16.39,16.41,16.43,16.45,16.47,16.49,16.51,16.53,16.55,16.57,16.59,16.61,16.63,16.65,16.68,16.70,16.72,16.74,16.76,16.78,16.80,16.83,16.84,16.87,16.89,16.91,16.93,16.95,16.97,17.00,17.02,17.04,17.07,17.09,17.12,17.13,17.16,17.18,17.20,17.22,17.24,17.26,17.29,17.31,17.33,17.35,17.38,17.40,17.42,17.45,17.47,17.50,17.52,17.54,17.57,17.59,17.62,17.64,17.67,17.70,17.72,17.75,17.77,17.80,17.82,17.85,17.88,17.90,17.92,17.95,17.97,18.00,18.02,18.05,18.08,18.10,18.13,18.15,18.18,18.21,18.23,18.26,18.29,18.32,18.34,18.37,18.40,18.42,18.45,18.48,18.51,18.53,18.56,18.59,18.62,18.65,18.67,18.70,18.73,18.76,18.79,18.82,18.85,18.88,18.91,18.93,18.97,19.00,19.02,19.06,19.09,19.11,19.14,19.17,19.21,19.23,19.26,19.30,19.32,19.35,19.39,19.41,19.45,19.48,19.51,19.54,19.58,19.61,19.64,19.67,19.71,19.74,19.77,19.81,19.84,19.87,19.91,19.94,19.97,20.01,20.04,20.07,20.11,20.14,20.18,20.21,20.24,20.28,20.32,20.35,20.39,20.42,20.46,20.50,20.53,20.57,20.61,20.64,20.67,20.71,20.75,20.78,20.81,20.85,20.89,20.93,20.96,21.00,21.04,21.08,21.11,21.16,21.20,21.23,21.27,21.31,21.35,21.40,21.43,21.47,21.51,21.55,21.59,21.63,21.68,21.71,21.75,21.79,21.83,21.88,21.92,21.96,22.00,22.05,22.09,22.13,22.18,22.22,22.26,22.30,22.35,22.39,22.43,22.48,22.52,22.57,22.61,22.66,22.71,22.76,22.80,22.85,22.90,22.95,22.99,23.04,23.09,23.14,23.19,23.23,23.28,23.33,23.38,23.42,23.47,23.52,23.57,23.61,23.66,23.71,23.76,23.80,23.85,23.90,23.96,24.01,24.06,24.11,24.17,24.22,24.27,24.32,24.38,24.43,24.49,24.54,24.59,24.64,24.70,24.76,24.81,24.87,24.92,24.98,25.04,25.09,25.15,25.21,25.27,25.32,25.38,25.44,25.50,25.56,25.61,25.68,25.74,25.80,25.86,25.92,25.98,26.04,26.10,26.17,26.23,26.29,26.35,26.41,26.48,26.54,26.61,26.68,26.74,26.80,26.87,26.94,27.00,27.07,27.14,27.21,27.27,27.34,27.41,27.48,27.55,27.62,27.69,27.76,27.83,27.90,27.97,28.05,28.12,28.19,28.27,28.34,28.42,28.49,28.57,28.65,28.72,28.80,28.87,28.95,29.03,29.11,29.18,29.26,29.34,29.42,29.50,29.58,29.66,29.75,29.83,29.91,30.00,30.08,30.16,30.25,30.34,30.42,30.51,30.60,30.69,30.77,30.86,30.95,31.04,31.13,31.22,31.31,31.40,31.49,31.59,31.68,31.78,31.87,31.97,32.06,32.16,32.26,32.36,32.46,32.56,32.66,32.76,32.86,32.97,33.07,33.18,33.28,33.39,33.49,33.60,33.71,33.82,33.92,34.03,34.14,34.25,34.36,34.47,34.59,34.70,34.82,34.93,35.05,35.16,35.28,35.40,35.52,35.64,35.76,35.88,36.01,36.13,36.25,36.38,36.51,36.63,36.76,36.89,37.02,37.16,37.29,37.42,37.55,37.69,37.82,37.96,38.10,38.24,38.38,38.52,38.66,38.80,38.95,39.10,39.24,39.39,39.54,39.69,39.84,40.00,40.15,40.31,40.47,40.62,40.78,40.94,41.10,41.27,41.43,41.60,41.77,41.94,42.10,42.28,42.45,42.63,42.80,42.98,43.16,43.34,43.52,43.70,43.89,44.07,44.26,44.45,44.64,44.84,45.03,45.23,45.43,45.63,45.84,46.04,46.25,46.46,46.67,46.88,47.09,47.31,47.53,47.75,47.97,48.19,48.42,48.65,48.87,49.11,49.35,49.58,49.82,50.06,50.31,50.56,50.81,51.06,51.32,51.58,51.84,52.10,52.36,52.63,52.90,53.18,53.45,53.74,54.02,54.31,54.60,54.89,55.19,55.49,55.79,56.10,56.41,56.72,57.03,57.35,57.68,58.00,58.33,58.67,59.01,59.35,59.70,60.05,60.40,60.77,61.13,61.50,61.88,62.26,62.64,63.03,63.42,63.81,64.21,64.62,65.03,65.45,65.87,66.30,66.73,67.17,67.62,68.07,68.53,68.99,69.46,69.94,70.43,70.92,71.42,71.92,72.44,72.95,73.48,74.02,74.56,75.11,75.67,76.23,76.81,77.40,77.99,78.60,79.21,79.83,80.47,81.11,81.77,82.44,83.11,83.80,84.50,85.22,85.94,86.68,87.43,88.19,88.96,89.75,90.56,91.37,92.20,93.05,93.91,94.79,95.69,96.52,96.95,97.39,97.84,98.29];

    ind = find(min(abs(nA_list-fund)) == abs(nA_list-fund),1); 

    ind = 2049 - ind; 
    
    M1 = M1_list(ind);
    M2 = M2_list(ind);
    M3 = M3_list(ind);  
    
end
% --------------------------------------------------------------

% --------------------------------------------------------------
function [M1,M2,M3,nVs,duty] = b12toMod(bit)

    M1_list = [250,249,249,249,248,248,248,247,247,247,246,246,246,246,245,245,245,244,244,244,243,243,243,242,242,242,241,241,241,240,240,240,239,239,239,239,238,238,238,237,237,237,236,236,236,235,235,235,234,234,234,233,233,233,232,232,232,232,231,231,231,230,230,230,229,229,229,228,228,228,227,227,227,226,226,226,225,225,225,224,224,224,224,223,223,223,222,222,222,221,221,221,220,220,220,219,219,219,218,218,218,217,217,217,217,216,216,216,215,215,215,214,214,214,213,213,213,212,212,212,211,211,211,210,210,210,209,209,209,209,208,208,208,207,207,207,206,206,206,205,205,205,204,204,204,203,203,203,202,202,202,202,201,201,201,200,200,200,199,199,199,198,198,198,197,197,197,196,196,196,195,195,195,195,194,194,194,193,193,193,192,192,192,191,191,191,190,190,190,189,189,189,188,188,188,187,187,187,187,186,186,186,185,185,185,184,184,184,183,183,183,182,182,182,181,181,181,180,180,180,180,179,179,179,178,178,178,177,177,177,176,176,176,175,175,175,174,174,174,173,173,173,173,172,172,172,171,171,171,170,170,170,169,169,169,168,168,168,167,167,167,166,166,166,165,165,165,165,164,164,164,163,163,163,162,162,162,161,161,161,160,160,160,159,159,159,158,158,158,158,157,157,157,156,156,156,155,155,155,154,154,154,153,153,153,152,152,152,151,151,151,151,150,150,150,149,149,149,148,148,148,147,147,147,146,146,146,145,145,145,144,144,144,143,143,143,143,142,142,142,141,141,141,140,140,140,139,139,139,138,138,138,137,137,137,136,136,136,136,135,135,135,134,134,134,133,133,133,132,132,132,131,131,131,130,130,130,129,129,129,128,128,128,128,127,127,127,126,126,126,125,125,125,124,124,124,123,123,123,122,122,122,122,121,121,121,120,120,119,119,119,118,118,118,117,117,116,116,116,115,115,115,114,114,114,113,113,113,112,112,112,111,111,110,110,110,109,109,109,108,108,107,107,107,106,106,106,105,105,105,104,104,104,103,103,102,102,102,101,101,101,100,100,100,99,99,99,98,98,97,97,97,96,96,96,95,95,95,94,94,94,93,93,92,92,92,91,91,91,90,90,90,89,89,89,88,88,87,87,87,86,86,86,85,85,85,84,84,84,83,83,82,82,82,81,81,81,80,80,80,79,79,78,78,78,77,77,77,76,76,76,75,75,75,74,74,73,73,73,72,72,72,71,71,71,70,70,70,69,69,68,68,68,67,67,67,66,66,66,65,65,65,64,64,64,63,63,62,62,62,61,61,61,60,60,60,59,59,58,58,58,57,57,57,56,56,56,55,55,55,55,55,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,48,48,48,48,48,48,48,48,48,48,48,48,48,48,48,48,48,48,48,48,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,44,44,44,44,44,44,44,44,44,44,44,44,44,44,44,44,44,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,41,41,41,41,41,41,41,41,41,41,41,41,41,41,41,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,38,38,38,38,38,38,38,38,38,38,38,38,38,38,38,37,37,37,37,37,37,37,37,37,37,37,37,37,37,36,36,36,36,36,36,36,36,36,36,36,36,36,35,35,35,35,35,35,35,35,35,35,35,35,35,35,34,34,34,34,34,34,34,34,34,34,34,34,34,34,33,33,33,33,33,33,33,33,33,33,33,33,33,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,4,4,4,4,4,4,4,4,4,4,4,4,4,4];
    M2_list = [250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,249,249,249,248,248,248,247,247,247,247,246,246,246,245,245,245,244,244,244,243,243,243,242,242,242,241,241,241,241,240,240,240,239,239,239,238,238,238,237,237,237,236,236,236,236,235,235,235,234,234,234,233,233,233,232,232,232,231,231,231,230,230,230,230,229,229,229,228,228,228,227,227,227,226,226,226,225,225,225,224,224,224,224,223,223,223,222,222,222,221,221,221,220,220,220,219,219,219,218,218,218,218,217,217,217,216,216,216,215,215,215,214,214,214,213,213,213,213,212,212,212,211,211,211,210,210,210,209,209,209,208,208,208,207,207,207,207,206,206,206,205,205,205,204,204,204,203,203,203,202,202,202,201,201,201,201,200,200,200,199,199,199,198,198,198,197,197,197,196,196,196,195,195,195,195,194,194,194,193,193,193,192,192,192,191,191,191,190,190,190,190,189,189,189,188,188,188,187,187,187,187,187,186,186,186,186,185,185,184,184,184,184,183,183,183,182,182,182,182,181,181,180,180,180,180,179,179,178,178,178,178,177,177,177,176,176,176,176,175,175,175,174,174,174,173,173,173,172,172,172,171,171,171,171,170,170,170,169,169,169,168,168,168,167,167,167,167,166,166,166,165,165,165,164,164,164,164,163,163,163,162,162,162,161,161,161,161,160,160,160,159,159,159,159,158,158,157,157,157,156,156,156,156,155,155,155,155,154,154,154,153,153,153,153,152,152,151,151,151,151,150,150,150,149,149,149,149,148,148,148,147,147,147,147,146,146,146,145,145,145,144,144,144,144,143,143,143,142,142,142,141,141,141,141,140,140,140,139,139,139,139,138,138,138,137,137,137,136,136,136,136,135,135,135,134,134,134,134,133,133,133,132,132,132,132,131,131,131,130,130,130,129,129,129,129,128,128,128,127,127,127,127,126,126,126,125,125,125,125,124,124,124,124,123,123,123,122,122,122,122,121,121,121,120,120,120,119,119,119,119,118,118,118,117,117,117,117,116,116,116,115,115,115,115,114,114,114,113,113,113,113,112,112,112,111,111,111,111,110,110,110,110,109,109,109,108,108,108,108,107,107,107,106,106,106,106,105,105,105,104,104,104,104,103,103,103,103,103,103,103,103,103,103,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,101,101,101,101,101,101,101,101,101,101,101,101,101,101,101,101,101,101,101,101,101,101,101,101,101,101,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,98,98,98,98,98,98,98,98,98,98,98,98,98,98,98,98,98,97,97,97,97,97,97,97,97,97,97,97,97,97,97,97,97,97,96,96,96,96,96,96,96,96,96,96,96,96,96,96,96,95,95,95,95,95,95,95,95,95,95,95,95,95,95,95,94,94,94,94,94,94,94,94,94,94,94,94,94,94,93,93,93,93,93,93,93,93,93,93,93,93,93,92,92,92,92,92,92,92,92,92,92,92,92,92,91,91,91,91,91,91,91,91,91,91,91,91,91,90,90,90,90,90,90,90,90,90,90,90,90,89,89,89,89,89,89,89,89,89,89,89,89,88,88,88,88,88,88,88,88,88,88,88,88,87,87,87,87,87,87,87,87,87,87,87,87,86,86,86,86,86,86,86,86,86,86,86,85,85,85,85,85,85,85,85,85,85,85,84,84,84,84,84,84,84,84,84,84,84,83,83,83,83,83,83,83,83,83,83,83,82,82,82,82,82,82,82,82,82,82,82,81,81,81,81,81,81,81,81,81,81,81,80,80,80,80,80,80,80,80,80,80,79,79,79,79,79,79,79,79,79,79,78,78,78,78,78,78,78,78,78,78,78,77,77,77,77,77,77,77,77,77,77,76,76,76,76,76,76,76,76,76,76,75,75,75,75,75,75,75,75,75,75,74,74,74,74,74,74,74,74,74,74,73,73,73,73,73,73,73,73,73,73,72,72,72,72,72,72,72,72,72,72,71,71,71,71,71,71,71,71,71,71,70,70,70,70,70,70,70,70,70,69,69,69,69,69,69,69,69,69,69,68,68,68,68,68,68,68,68,68,67,67,67,67,67,67,67,67,67,67,66,66,66,66,66,66,66,66,66,65,65,65,65,65,65,65,65,65,64,64,64,64,64,64,64,64,64,64,63,63,63,63,63,63,63,63,63,62,62,62,62,62,62,62,62,62,61,61,61,61,61,61,61,61,61,61,60,60,60,60,60,60,60,60,60,59,59,59,59,59,59,59,59,59,58,58,58,58,58,58,58,58,58,57,57,57,57,57,57,57,57,57,56,56,56,56,56,56,56,56,56,55,55,55,55,55,55,55,55,55,54,54,54,54,54,54,54,54,54,53,53,53,53,53,53,53,53,53,52,52,52,52,52,52,52,52,52,51,51,51,51,51,51,51,51,51,50,50,50,50,50,50,50,50,50,49,49,49,49,49,49,49,49,49,48,48,48,48,48,48,48,48,48,47,47,47,47,47,47,47,47,46,46,46,46,46,46,46,46,46,45,45,45,45,45,45,45,45,45,44,44,44,44,44,44,44,44,44,43,43,43,43,43,43,43,43,43,42,42,42,42,42,42,42,42,41,41,41,41,41,41,41,41,41,40,40,40,40,40,40,40,40,40,39,39,39,39,39,39,39,39,38,38,38,38,38,38,38,38,38,37,37,37,37,37,37,37,37,37,36,36,36,36,36,36,36,36,35,35,35,35,35,35,35,35,35,34,34,34,34,34,34,34,34,34,33,33,33,33,33,33,33,33,32,32,32,32,32,32,32,32,32,31,31,31,31,31,31,31,31,30,30,30,30,30,30,30,30,30,29,29,29,29,29,29,29,29,29,28,28,28,28,28,28,28,28,27,27,27,27,27,27,27,27,27,26,26,26,26,26,26,26,26,25,25,25,25,25,25,25,25,25,24,24,24,24,24,24,24,24,23,23,23,23,23,23,23,23,23,22,22,22,22,22,22,22,22,21,21,21,21,21,21,21,21,21,20,20,20,20,20,20,20,20,19,19,19,19,19,19,19,19,19,18,18,18,18,18,18,18,18,17,17,17,17,17,17,17,17,17,16,16,16,16,16,16,16,16,15,15,15,15,15,15,15,15,14,14,14,14,14,14,14,14,14,13,13,13,13,13,13,13,13,12,12,12,12,12,12];
    M3_list = [250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,249,249,249,249,248,248,248,248,247,247,247,246,246,246,246,245,245,245,245,244,244,244,243,243,243,243,242,242,242,241,241,241,241,240,240,240,240,239,239,239,238,238,238,238,237,237,237,237,236,236,236,235,235,235,235,234,234,234,234,233,233,233,232,232,232,232,231,231,231,230,230,230,230,229,229,229,229,228,228,228,227,227,227,227,226,226,226,226,225,225,225,224,224,224,224,223,223,223,223,222,222,222,221,221,221,221,220,220,220,219,219,219,219,218,218,218,218,217,217,217,216,216,216,216,215,215,215,215,214,214,214,213,213,213,213,212,212,212,212,211,211,211,210,210,210,210,209,209,209,208,208,208,208,207,207,207,207,206,206,206,205,205,205,205,204,204,204,204,203,203,203,202,202,202,202,201,201,201,201,200,200,200,199,199,199,199,198,198,198,197,197,197,197,196,196,196,196,195,195,195,195,195,195,194,194,194,194,193,193,193,192,192,192,192,191,191,191,191,190,190,190,189,189,189,189,188,188,188,187,187,187,187,186,186,186,186,185,185,185,185,184,184,184,184,183,183,183,182,182,182,182,182,181,181,181,180,180,180,180,179,179,179,179,178,178,178,178,177,177,177,177,176,176,176,176,175,175,175,175,174,174,174,174,173,173,173,173,172,172,172,172,171,171,171,171,170,170,170,170,169,169,169,169,168,168,168,168,167,167,167,167,166,166,166,166,165,165,165,165,165,164,164,164,163,163,163,163,163,162,162,162,162,161,161,161,161,160,160,160,160,159,159,159,159,158,158,158,158,157,157,157,157,156,156,156,156,156,155,155,155,155,154,154,154,154,153,153,153,153,152,152,152,152,152,151,151,151,151,150,150,150,150,149,149,149,149,148,148,148,148,148,147,147,147,147,146,146,146,146,145,145,145,145,144,144,144,144,144,143,143,143,143,142,142,142,142,142,141,141,141,141,140,140,140,140,139,139,139,139,139,138,138,138,138,137,137,137,137,136,136,136,136,136,135,135,135,135,134,134,134,134,134,133,133,133,133,132,132,132,132,132,131,131,131,131,130,130,130,130,129,129,129,129,129,128,128,128,128,128,127,127,127,127,126,126,126,126,126,125,125,125,125,124,124,124,124,123,123,123,123,123,122,122,122,122,122,121,121,121,121,120,120,120,120,120,119,119,119,119,118,118,118,118,118,117,117,117,117,117,116,116,116,116,115,115,115,115,115,114,114,114,114,113,113,113,113,113,112,112,112,112,112,111,111,111,111,110,110,110,110,110,109,109,109,109,108,108,108,108,108,107,107,107,107,107,106,106,106,106,105,105,105,105,105,104,104,104,104,104,103,103,103,103,103,102,102,102,102,101,101,101,101,101,100,100,100,100,100,99,99,99,99,98,98,98,98,98,97,97,97,97,97,96,96,96,96,96,95,95,95,95,94,94,94,94,94,93,93,93,93,93,92,92,92,92,91,91,91,91,91,90,90,90,90,90,89,89,89,89,89,88,88,88,88,88,87,87,87,87,87,86,86,86,86,85,85,85,85,85,84,84,84,84,84,83,83,83,83,83,82,82,82,82,82,81,81,81,81,81,80,80,80,80,79,79,79,79,79,78,78,78,78,78,77,77,77,77,77,76,76,76,76,76,75,75,75,75,75,74,74,74,74,74,73,73,73,73,72,72,72,72,72,71,71,71,71,71,70,70,70,70,70,69,69,69,69,69,68,68,68,68,68,67,67,67,67,67,66,66,66,66,66,65,65,65,65,65,64,64,64,64,64,63,63,63,63,63,62,62,62,62,61,61,61,61,61,60,60,60,60,60,59,59,59,59,59,58,58,58,58,58,57,57,57,57,57,56,56,56,56,56,55,55,55,55,55,54,54,54,54,54,53,53,53,53,53,52,52,52,52,52,51,51,51,51,51,50,50,50,50,50,49,49,49,49,49,48,48,48,48,48,47,47,47,47,47,46,46,46,46,46,45,45,45,45,45,44,44,44,44,44,43,43,43,43,43,42,42,42,42,42,41,41,41,41,41,40,40,40,40,40,39,39,39,39,39,38,38,38,38,38,37,37,37,37,37,36,36,36,36,36,35,35,35,35,35,34,34,34,34,34,33,33,33,33,33,32,32,32,32,32,31,31,31,31,31,30,30,30,30,30,29,29,29,29,29,28,28,28,28,28,27,27,27,27,27,26,26,26,26,26,25,25,25,25,25,24,24,24,24,24,23,23,23,23,23,22,22,22,22,22,21,21,21,21,21,21,20,20,20];VrecFund = [2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.51, 2.51, 2.52, 2.52, 2.53, 2.53, 2.54, 2.54, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.56, 2.56, 2.57, 2.57, 2.58, 2.58, 2.59, 2.59, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.61, 2.61, 2.62, 2.62, 2.63, 2.63, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.65, 2.65, 2.66, 2.66, 2.67, 2.67, 2.68, 2.68, 2.69, 2.69, 2.69, 2.69, 2.69, 2.69, 2.69, 2.69, 2.69, 2.69, 2.69, 2.70, 2.70, 2.70, 2.71, 2.71, 2.72, 2.72, 2.73, 2.73, 2.74, 2.74, 2.74, 2.74, 2.74, 2.74, 2.74, 2.74, 2.74, 2.75, 2.75, 2.76, 2.76, 2.77, 2.77, 2.78, 2.78, 2.79, 2.79, 2.79, 2.79, 2.79, 2.79, 2.79, 2.80, 2.80, 2.80, 2.81, 2.81, 2.82, 2.82, 2.83, 2.83, 2.83, 2.83, 2.83, 2.84, 2.84, 2.85, 2.85, 2.86, 2.86, 2.87, 2.87, 2.88, 2.88, 2.88, 2.88, 2.88, 2.89, 2.89, 2.90, 2.90, 2.90, 2.91, 2.91, 2.92, 2.92, 2.93, 2.93, 2.93, 2.93, 2.94, 2.94, 2.95, 2.95, 2.96, 2.96, 2.97, 2.97, 2.98, 2.98, 2.98, 2.99, 2.99, 3.00, 3.00, 3.00, 3.01, 3.01, 3.02, 3.02, 3.03, 3.03, 3.04, 3.04, 3.05, 3.05, 3.06, 3.06, 3.07, 3.08, 3.08, 3.09, 3.09, 3.10, 3.10, 3.10, 3.11, 3.11, 3.12, 3.13, 3.13, 3.14, 3.14, 3.15, 3.15, 3.16, 3.16, 3.17, 3.18, 3.18, 3.19, 3.19, 3.20, 3.20, 3.20, 3.21, 3.22, 3.23, 3.23, 3.24, 3.24, 3.25, 3.26, 3.27, 3.28, 3.28, 3.29, 3.29, 3.30, 3.30, 3.31, 3.32, 3.32, 3.33, 3.33, 3.34, 3.34, 3.35, 3.36, 3.37, 3.38, 3.38, 3.39, 3.40, 3.40, 3.41, 3.42, 3.43, 3.43, 3.44, 3.45, 3.46, 3.47, 3.47, 3.48, 3.48, 3.49, 3.50, 3.51, 3.52, 3.52, 3.53, 3.54, 3.55, 3.56, 3.57, 3.57, 3.58, 3.59, 3.59, 3.60, 3.61, 3.62, 3.63, 3.64, 3.65, 3.66, 3.67, 3.68, 3.69, 3.69, 3.70, 3.71, 3.72, 3.73, 3.74, 3.75, 3.76, 3.77, 3.78, 3.79, 3.80, 3.81, 3.82, 3.83, 3.84, 3.85, 3.86, 3.87, 3.88, 3.89, 3.89, 3.91, 3.92, 3.93, 3.94, 3.95, 3.97, 3.98, 3.99, 4.00, 4.01, 4.02, 4.03, 4.05, 4.06, 4.07, 4.08, 4.09, 4.11, 4.12, 4.13, 4.14, 4.16, 4.17, 4.18, 4.19, 4.20, 4.21, 4.23, 4.24, 4.26, 4.27, 4.28, 4.29, 4.31, 4.32, 4.34, 4.35, 4.36, 4.38, 4.39, 4.40, 4.42, 4.43, 4.45, 4.46, 4.48, 4.49, 4.51, 4.52, 4.54, 4.56, 4.57, 4.59, 4.60, 4.62, 4.63, 4.65, 4.67, 4.68, 4.70, 4.71, 4.73, 4.75, 4.77, 4.79, 4.80, 4.82, 4.84, 4.86, 4.88, 4.89, 4.91, 4.93, 4.95, 4.97, 4.99, 5.01, 5.03, 5.05, 5.07, 5.09, 5.11, 5.13, 5.15, 5.17, 5.19, 5.21, 5.24, 5.26, 5.28, 5.30, 5.32, 5.35, 5.37, 5.39, 5.41, 5.43, 5.46, 5.48, 5.51, 5.53, 5.56, 5.58, 5.61, 5.64, 5.67, 5.69, 5.71, 5.74, 5.77, 5.79, 5.82, 5.85, 5.88, 5.90, 5.93, 5.97, 5.99, 6.02, 6.05, 6.08, 6.11, 6.14, 6.17, 6.20, 6.24, 6.27, 6.30, 6.33, 6.37, 6.40, 6.43, 6.47, 6.51, 6.54, 6.58, 6.61, 6.65, 6.69, 6.73, 6.77, 6.80, 6.84, 6.88, 6.92, 6.97, 7.00, 7.04, 7.08, 7.12, 7.17, 7.21, 7.25, 7.29, 7.32, 7.35, 7.38, 7.40, 7.43, 7.45, 7.46, 7.47, 7.47, 7.48, 7.48, 7.49, 7.49, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.51, 7.51, 7.52, 7.52, 7.53, 7.53, 7.54, 7.54, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.56, 7.56, 7.57, 7.57, 7.57, 7.58, 7.58, 7.59, 7.59, 7.59, 7.59, 7.59, 7.59, 7.59, 7.59, 7.59, 7.60, 7.60, 7.61, 7.61, 7.62, 7.62, 7.63, 7.63, 7.64, 7.64, 7.64, 7.64, 7.64, 7.64, 7.64, 7.65, 7.65, 7.66, 7.66, 7.67, 7.67, 7.67, 7.68, 7.68, 7.69, 7.69, 7.69, 7.69, 7.70, 7.70, 7.71, 7.71, 7.72, 7.72, 7.73, 7.73, 7.74, 7.74, 7.74, 7.74, 7.75, 7.75, 7.76, 7.76, 7.77, 7.77, 7.77, 7.78, 7.78, 7.79, 7.79, 7.80, 7.80, 7.81, 7.81, 7.82, 7.82, 7.83, 7.83, 7.84, 7.84, 7.85, 7.85, 7.86, 7.86, 7.87, 7.87, 7.87, 7.88, 7.89, 7.89, 7.90, 7.90, 7.91, 7.91, 7.92, 7.92, 7.93, 7.94, 7.94, 7.95, 7.95, 7.96, 7.96, 7.97, 7.97, 7.98, 7.99, 7.99, 8.00, 8.00, 8.01, 8.01, 8.02, 8.03, 8.04, 8.04, 8.05, 8.05, 8.06, 8.07, 8.07, 8.08, 8.09, 8.09, 8.10, 8.10, 8.11, 8.12, 8.13, 8.14, 8.14, 8.15, 8.15, 8.16, 8.17, 8.18, 8.18, 8.19, 8.19, 8.20, 8.21, 8.22, 8.23, 8.23, 8.24, 8.24, 8.25, 8.26, 8.27, 8.28, 8.28, 8.29, 8.30, 8.31, 8.32, 8.33, 8.33, 8.34, 8.35, 8.36, 8.37, 8.37, 8.38, 8.38, 8.39, 8.40, 8.41, 8.42, 8.43, 8.44, 8.45, 8.46, 8.47, 8.47, 8.48, 8.49, 8.50, 8.51, 8.52, 8.52, 8.53, 8.54, 8.55, 8.56, 8.57, 8.58, 8.59, 8.60, 8.61, 8.62, 8.63, 8.64, 8.65, 8.66, 8.67, 8.67, 8.68, 8.69, 8.70, 8.71, 8.72, 8.73, 8.74, 8.75, 8.76, 8.77, 8.78, 8.79, 8.80, 8.81, 8.82, 8.83, 8.84, 8.85, 8.87, 8.87, 8.88, 8.89, 8.91, 8.92, 8.93, 8.94, 8.95, 8.97, 8.97, 8.98, 9.00, 9.01, 9.02, 9.03, 9.04, 9.06, 9.07, 9.07, 9.08, 9.10, 9.11, 9.12, 9.13, 9.14, 9.16, 9.17, 9.18, 9.19, 9.21, 9.22, 9.23, 9.25, 9.26, 9.27, 9.28, 9.30, 9.31, 9.32, 9.33, 9.35, 9.36, 9.37, 9.38, 9.40, 9.41, 9.42, 9.44, 9.45, 9.46, 9.48, 9.49, 9.50, 9.52, 9.53, 9.55, 9.56, 9.57, 9.59, 9.60, 9.62, 9.63, 9.65, 9.66, 9.67, 9.69, 9.70, 9.72, 9.73, 9.75, 9.76, 9.77, 9.79, 9.81, 9.82, 9.84, 9.86, 9.87, 9.89, 9.90, 9.92, 9.94, 9.95, 9.96, 9.98,10.00,10.01,10.03,10.05,10.06,10.08,10.09,10.11,10.13,10.14,10.16,10.18,10.19,10.21,10.23,10.25,10.26,10.28,10.30,10.32,10.34,10.36,10.37,10.39,10.41,10.43,10.45,10.46,10.48,10.50,10.52,10.54,10.56,10.57,10.59,10.61,10.63,10.65,10.67,10.69,10.71,10.73,10.75,10.77,10.79,10.81,10.83,10.85,10.87,10.89,10.91,10.93,10.96,10.97,10.99,11.02,11.04,11.06,11.08,11.10,11.12,11.15,11.16,11.19,11.21,11.24,11.26,11.28,11.31,11.33,11.35,11.38,11.40,11.43,11.45,11.47,11.50,11.52,11.55,11.57,11.59,11.62,11.64,11.66,11.69,11.71,11.74,11.76,11.78,11.81,11.83,11.85,11.88,11.91,11.94,11.96,11.99,12.01,12.04,12.07,12.10,12.13,12.15,12.18,12.21,12.24,12.26,12.29,12.32,12.35,12.37,12.39,12.41,12.43,12.44,12.45,12.46,12.47,12.48,12.48,12.49,12.49,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.51,12.51,12.52,12.52,12.53,12.53,12.54,12.54,12.55,12.55,12.55,12.55,12.55,12.55,12.55,12.55,12.55,12.55,12.55,12.55,12.55,12.56,12.56,12.57,12.57,12.58,12.58,12.59,12.59,12.59,12.59,12.59,12.59,12.59,12.60,12.60,12.61,12.61,12.62,12.62,12.63,12.63,12.64,12.64,12.64,12.64,12.65,12.65,12.65,12.66,12.66,12.67,12.67,12.68,12.68,12.69,12.69,12.69,12.70,12.70,12.71,12.71,12.72,12.72,12.73,12.73,12.74,12.74,12.75,12.75,12.75,12.76,12.76,12.77,12.77,12.78,12.79,12.79,12.80,12.80,12.81,12.81,12.82,12.82,12.83,12.84,12.85,12.85,12.85,12.86,12.86,12.87,12.88,12.89,12.89,12.90,12.90,12.91,12.91,12.92,12.93,12.94,12.94,12.95,12.95,12.95,12.96,12.97,12.98,12.99,12.99,13.00,13.00,13.01,13.02,13.03,13.04,13.04,13.05,13.05,13.05,13.06,13.07,13.08,13.09,13.09,13.10,13.10,13.11,13.12,13.13,13.14,13.14,13.15,13.15,13.16,13.17,13.18,13.19,13.19,13.20,13.21,13.22,13.23,13.24,13.24,13.25,13.26,13.27,13.28,13.28,13.29,13.30,13.31,13.32,13.33,13.33,13.34,13.35,13.36,13.37,13.38,13.38,13.39,13.40,13.41,13.42,13.43,13.44,13.44,13.45,13.46,13.47,13.48,13.49,13.50,13.51,13.52,13.53,13.54,13.55,13.56,13.57,13.58,13.59,13.60,13.61,13.62,13.63,13.64,13.64,13.65,13.66,13.67,13.68,13.69,13.70,13.72,13.73,13.74,13.74,13.75,13.77,13.78,13.79,13.80,13.81,13.82,13.83,13.84,13.84,13.86,13.87,13.88,13.89,13.90,13.92,13.93,13.94,13.94,13.96,13.97,13.98,13.99,14.00,14.01,14.02,14.03,14.04,14.06,14.07,14.08,14.09,14.11,14.12,14.13,14.14,14.15,14.16,14.17,14.19,14.20,14.21,14.22,14.24,14.25,14.26,14.27,14.28,14.30,14.31,14.32,14.34,14.35,14.36,14.37,14.39,14.40,14.41,14.43,14.44,14.45,14.46,14.47,14.49,14.50,14.51,14.53,14.54,14.55,14.56,14.58,14.59,14.61,14.62,14.64,14.64,14.66,14.67,14.69,14.70,14.72,14.73,14.74,14.76,14.77,14.79,14.80,14.82,14.83,14.84,14.86,14.87,14.89,14.90,14.92,14.94,14.95,14.96,14.98,14.99,15.01,15.02,15.04,15.05,15.06,15.08,15.09,15.11,15.13,15.14,15.15,15.17,15.19,15.20,15.22,15.23,15.25,15.27,15.28,15.30,15.32,15.33,15.34,15.36,15.38,15.39,15.41,15.43,15.44,15.46,15.48,15.50,15.52,15.53,15.55,15.57,15.58,15.60,15.62,15.63,15.65,15.66,15.68,15.70,15.72,15.73,15.75,15.77,15.79,15.81,15.83,15.84,15.86,15.88,15.90,15.91,15.93,15.95,15.97,15.99,16.00,16.02,16.04,16.06,16.08,16.10,16.12,16.13,16.15,16.17,16.19,16.21,16.23,16.25,16.27,16.29,16.31,16.33,16.35,16.37,16.39,16.41,16.43,16.45,16.47,16.49,16.51,16.53,16.55,16.57,16.59,16.61,16.63,16.65,16.68,16.70,16.72,16.74,16.76,16.78,16.80,16.83,16.84,16.87,16.89,16.91,16.93,16.95,16.97,17.00,17.02,17.04,17.07,17.09,17.12,17.13,17.16,17.18,17.20,17.22,17.24,17.26,17.29,17.31,17.33,17.35,17.38,17.40,17.42,17.45,17.47,17.50,17.52,17.54,17.57,17.59,17.62,17.64,17.67,17.70,17.72,17.75,17.77,17.80,17.82,17.85,17.88,17.90,17.92,17.95,17.97,18.00,18.02,18.05,18.08,18.10,18.13,18.15,18.18,18.21,18.23,18.26,18.29,18.32,18.34,18.37,18.40,18.42,18.45,18.48,18.51,18.53,18.56,18.59,18.62,18.65,18.67,18.70,18.73,18.76,18.79,18.82,18.85,18.88,18.91,18.93,18.97,19.00,19.02,19.06,19.09,19.11,19.14,19.17,19.21,19.23,19.26,19.30,19.32,19.35,19.39,19.41,19.45,19.48,19.51,19.54,19.58,19.61,19.64,19.67,19.71,19.74,19.77,19.81,19.84,19.87,19.91,19.94,19.97,20.01,20.04,20.07,20.11,20.14,20.18,20.21,20.24,20.28,20.32,20.35,20.39,20.42,20.46,20.50,20.53,20.57,20.61,20.64,20.67,20.71,20.75,20.78,20.81,20.85,20.89,20.93,20.96,21.00,21.04,21.08,21.11,21.16,21.20,21.23,21.27,21.31,21.35,21.40,21.43,21.47,21.51,21.55,21.59,21.63,21.68,21.71,21.75,21.79,21.83,21.88,21.92,21.96,22.00,22.05,22.09,22.13,22.18,22.22,22.26,22.30,22.35,22.39,22.43,22.48,22.52,22.57,22.61,22.66,22.71,22.76,22.80,22.85,22.90,22.95,22.99,23.04,23.09,23.14,23.19,23.23,23.28,23.33,23.38,23.42,23.47,23.52,23.57,23.61,23.66,23.71,23.76,23.80,23.85,23.90,23.96,24.01,24.06,24.11,24.17,24.22,24.27,24.32,24.38,24.43,24.49,24.54,24.59,24.64,24.70,24.76,24.81,24.87,24.92,24.98,25.04,25.09,25.15,25.21,25.27,25.32,25.38,25.44,25.50,25.56,25.61,25.68,25.74,25.80,25.86,25.92,25.98,26.04,26.10,26.17,26.23,26.29,26.35,26.41,26.48,26.54,26.61,26.68,26.74,26.80,26.87,26.94,27.00,27.07,27.14,27.21,27.27,27.34,27.41,27.48,27.55,27.62,27.69,27.76,27.83,27.90,27.97,28.05,28.12,28.19,28.27,28.34,28.42,28.49,28.57,28.65,28.72,28.80,28.87,28.95,29.03,29.11,29.18,29.26,29.34,29.42,29.50,29.58,29.66,29.75,29.83,29.91,30.00,30.08,30.16,30.25,30.34,30.42,30.51,30.60,30.69,30.77,30.86,30.95,31.04,31.13,31.22,31.31,31.40,31.49,31.59,31.68,31.78,31.87,31.97,32.06,32.16,32.26,32.36,32.46,32.56,32.66,32.76,32.86,32.97,33.07,33.18,33.28,33.39,33.49,33.60,33.71,33.82,33.92,34.03,34.14,34.25,34.36,34.47,34.59,34.70,34.82,34.93,35.05,35.16,35.28,35.40,35.52,35.64,35.76,35.88,36.01,36.13,36.25,36.38,36.51,36.63,36.76,36.89,37.02,37.16,37.29,37.42,37.55,37.69,37.82,37.96,38.10,38.24,38.38,38.52,38.66,38.80,38.95,39.10,39.24,39.39,39.54,39.69,39.84,40.00,40.15,40.31,40.47,40.62,40.78,40.94,41.10,41.27,41.43,41.60,41.77,41.94,42.10,42.28,42.45,42.63,42.80,42.98,43.16,43.34,43.52,43.70,43.89,44.07,44.26,44.45,44.64,44.84,45.03,45.23,45.43,45.63,45.84,46.04,46.25,46.46,46.67,46.88,47.09,47.31,47.53,47.75,47.97,48.19,48.42,48.65,48.87,49.11,49.35,49.58,49.82,50.06,50.31,50.56,50.81,51.06,51.32,51.58,51.84,52.10,52.36,52.63,52.90,53.18,53.45,53.74,54.02,54.31,54.60,54.89,55.19,55.49,55.79,56.10,56.41,56.72,57.03,57.35,57.68,58.00,58.33,58.67,59.01,59.35,59.70,60.05,60.40,60.77,61.13,61.50,61.88,62.26,62.64,63.03,63.42,63.81,64.21,64.62,65.03,65.45,65.87,66.30,66.73,67.17,67.62,68.07,68.53,68.99,69.46,69.94,70.43,70.92,71.42,71.92,72.44,72.95,73.48,74.02,74.56,75.11,75.67,76.23,76.81,77.40,77.99,78.60,79.21,79.83,80.47,81.11,81.77,82.44,83.11,83.80,84.50,85.22,85.94,86.68,87.43,88.19,88.96,89.75,90.56,91.37,92.20,93.05,93.91,94.79,95.69,96.52,96.95,97.39,97.84,98.29];
    VrecFund = [2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.51, 2.51, 2.52, 2.52, 2.53, 2.53, 2.54, 2.54, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.56, 2.56, 2.57, 2.57, 2.58, 2.58, 2.59, 2.59, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.61, 2.61, 2.62, 2.62, 2.63, 2.63, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.65, 2.65, 2.66, 2.66, 2.67, 2.67, 2.68, 2.68, 2.69, 2.69, 2.69, 2.69, 2.69, 2.69, 2.69, 2.69, 2.69, 2.69, 2.69, 2.70, 2.70, 2.70, 2.71, 2.71, 2.72, 2.72, 2.73, 2.73, 2.74, 2.74, 2.74, 2.74, 2.74, 2.74, 2.74, 2.74, 2.74, 2.75, 2.75, 2.76, 2.76, 2.77, 2.77, 2.78, 2.78, 2.79, 2.79, 2.79, 2.79, 2.79, 2.79, 2.79, 2.80, 2.80, 2.80, 2.81, 2.81, 2.82, 2.82, 2.83, 2.83, 2.83, 2.83, 2.83, 2.84, 2.84, 2.85, 2.85, 2.86, 2.86, 2.87, 2.87, 2.88, 2.88, 2.88, 2.88, 2.88, 2.89, 2.89, 2.90, 2.90, 2.90, 2.91, 2.91, 2.92, 2.92, 2.93, 2.93, 2.93, 2.93, 2.94, 2.94, 2.95, 2.95, 2.96, 2.96, 2.97, 2.97, 2.98, 2.98, 2.98, 2.99, 2.99, 3.00, 3.00, 3.00, 3.01, 3.01, 3.02, 3.02, 3.03, 3.03, 3.04, 3.04, 3.05, 3.05, 3.06, 3.06, 3.07, 3.08, 3.08, 3.09, 3.09, 3.10, 3.10, 3.10, 3.11, 3.11, 3.12, 3.13, 3.13, 3.14, 3.14, 3.15, 3.15, 3.16, 3.16, 3.17, 3.18, 3.18, 3.19, 3.19, 3.20, 3.20, 3.20, 3.21, 3.22, 3.23, 3.23, 3.24, 3.24, 3.25, 3.26, 3.27, 3.28, 3.28, 3.29, 3.29, 3.30, 3.30, 3.31, 3.32, 3.32, 3.33, 3.33, 3.34, 3.34, 3.35, 3.36, 3.37, 3.38, 3.38, 3.39, 3.40, 3.40, 3.41, 3.42, 3.43, 3.43, 3.44, 3.45, 3.46, 3.47, 3.47, 3.48, 3.48, 3.49, 3.50, 3.51, 3.52, 3.52, 3.53, 3.54, 3.55, 3.56, 3.57, 3.57, 3.58, 3.59, 3.59, 3.60, 3.61, 3.62, 3.63, 3.64, 3.65, 3.66, 3.67, 3.68, 3.69, 3.69, 3.70, 3.71, 3.72, 3.73, 3.74, 3.75, 3.76, 3.77, 3.78, 3.79, 3.80, 3.81, 3.82, 3.83, 3.84, 3.85, 3.86, 3.87, 3.88, 3.89, 3.89, 3.91, 3.92, 3.93, 3.94, 3.95, 3.97, 3.98, 3.99, 4.00, 4.01, 4.02, 4.03, 4.05, 4.06, 4.07, 4.08, 4.09, 4.11, 4.12, 4.13, 4.14, 4.16, 4.17, 4.18, 4.19, 4.20, 4.21, 4.23, 4.24, 4.26, 4.27, 4.28, 4.29, 4.31, 4.32, 4.34, 4.35, 4.36, 4.38, 4.39, 4.40, 4.42, 4.43, 4.45, 4.46, 4.48, 4.49, 4.51, 4.52, 4.54, 4.56, 4.57, 4.59, 4.60, 4.62, 4.63, 4.65, 4.67, 4.68, 4.70, 4.71, 4.73, 4.75, 4.77, 4.79, 4.80, 4.82, 4.84, 4.86, 4.88, 4.89, 4.91, 4.93, 4.95, 4.97, 4.99, 5.01, 5.03, 5.05, 5.07, 5.09, 5.11, 5.13, 5.15, 5.17, 5.19, 5.21, 5.24, 5.26, 5.28, 5.30, 5.32, 5.35, 5.37, 5.39, 5.41, 5.43, 5.46, 5.48, 5.51, 5.53, 5.56, 5.58, 5.61, 5.64, 5.67, 5.69, 5.71, 5.74, 5.77, 5.79, 5.82, 5.85, 5.88, 5.90, 5.93, 5.97, 5.99, 6.02, 6.05, 6.08, 6.11, 6.14, 6.17, 6.20, 6.24, 6.27, 6.30, 6.33, 6.37, 6.40, 6.43, 6.47, 6.51, 6.54, 6.58, 6.61, 6.65, 6.69, 6.73, 6.77, 6.80, 6.84, 6.88, 6.92, 6.97, 7.00, 7.04, 7.08, 7.12, 7.17, 7.21, 7.25, 7.29, 7.32, 7.35, 7.38, 7.40, 7.43, 7.45, 7.46, 7.47, 7.47, 7.48, 7.48, 7.49, 7.49, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.51, 7.51, 7.52, 7.52, 7.53, 7.53, 7.54, 7.54, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.56, 7.56, 7.57, 7.57, 7.57, 7.58, 7.58, 7.59, 7.59, 7.59, 7.59, 7.59, 7.59, 7.59, 7.59, 7.59, 7.60, 7.60, 7.61, 7.61, 7.62, 7.62, 7.63, 7.63, 7.64, 7.64, 7.64, 7.64, 7.64, 7.64, 7.64, 7.65, 7.65, 7.66, 7.66, 7.67, 7.67, 7.67, 7.68, 7.68, 7.69, 7.69, 7.69, 7.69, 7.70, 7.70, 7.71, 7.71, 7.72, 7.72, 7.73, 7.73, 7.74, 7.74, 7.74, 7.74, 7.75, 7.75, 7.76, 7.76, 7.77, 7.77, 7.77, 7.78, 7.78, 7.79, 7.79, 7.80, 7.80, 7.81, 7.81, 7.82, 7.82, 7.83, 7.83, 7.84, 7.84, 7.85, 7.85, 7.86, 7.86, 7.87, 7.87, 7.87, 7.88, 7.89, 7.89, 7.90, 7.90, 7.91, 7.91, 7.92, 7.92, 7.93, 7.94, 7.94, 7.95, 7.95, 7.96, 7.96, 7.97, 7.97, 7.98, 7.99, 7.99, 8.00, 8.00, 8.01, 8.01, 8.02, 8.03, 8.04, 8.04, 8.05, 8.05, 8.06, 8.07, 8.07, 8.08, 8.09, 8.09, 8.10, 8.10, 8.11, 8.12, 8.13, 8.14, 8.14, 8.15, 8.15, 8.16, 8.17, 8.18, 8.18, 8.19, 8.19, 8.20, 8.21, 8.22, 8.23, 8.23, 8.24, 8.24, 8.25, 8.26, 8.27, 8.28, 8.28, 8.29, 8.30, 8.31, 8.32, 8.33, 8.33, 8.34, 8.35, 8.36, 8.37, 8.37, 8.38, 8.38, 8.39, 8.40, 8.41, 8.42, 8.43, 8.44, 8.45, 8.46, 8.47, 8.47, 8.48, 8.49, 8.50, 8.51, 8.52, 8.52, 8.53, 8.54, 8.55, 8.56, 8.57, 8.58, 8.59, 8.60, 8.61, 8.62, 8.63, 8.64, 8.65, 8.66, 8.67, 8.67, 8.68, 8.69, 8.70, 8.71, 8.72, 8.73, 8.74, 8.75, 8.76, 8.77, 8.78, 8.79, 8.80, 8.81, 8.82, 8.83, 8.84, 8.85, 8.87, 8.87, 8.88, 8.89, 8.91, 8.92, 8.93, 8.94, 8.95, 8.97, 8.97, 8.98, 9.00, 9.01, 9.02, 9.03, 9.04, 9.06, 9.07, 9.07, 9.08, 9.10, 9.11, 9.12, 9.13, 9.14, 9.16, 9.17, 9.18, 9.19, 9.21, 9.22, 9.23, 9.25, 9.26, 9.27, 9.28, 9.30, 9.31, 9.32, 9.33, 9.35, 9.36, 9.37, 9.38, 9.40, 9.41, 9.42, 9.44, 9.45, 9.46, 9.48, 9.49, 9.50, 9.52, 9.53, 9.55, 9.56, 9.57, 9.59, 9.60, 9.62, 9.63, 9.65, 9.66, 9.67, 9.69, 9.70, 9.72, 9.73, 9.75, 9.76, 9.77, 9.79, 9.81, 9.82, 9.84, 9.86, 9.87, 9.89, 9.90, 9.92, 9.94, 9.95, 9.96, 9.98,10.00,10.01,10.03,10.05,10.06,10.08,10.09,10.11,10.13,10.14,10.16,10.18,10.19,10.21,10.23,10.25,10.26,10.28,10.30,10.32,10.34,10.36,10.37,10.39,10.41,10.43,10.45,10.46,10.48,10.50,10.52,10.54,10.56,10.57,10.59,10.61,10.63,10.65,10.67,10.69,10.71,10.73,10.75,10.77,10.79,10.81,10.83,10.85,10.87,10.89,10.91,10.93,10.96,10.97,10.99,11.02,11.04,11.06,11.08,11.10,11.12,11.15,11.16,11.19,11.21,11.24,11.26,11.28,11.31,11.33,11.35,11.38,11.40,11.43,11.45,11.47,11.50,11.52,11.55,11.57,11.59,11.62,11.64,11.66,11.69,11.71,11.74,11.76,11.78,11.81,11.83,11.85,11.88,11.91,11.94,11.96,11.99,12.01,12.04,12.07,12.10,12.13,12.15,12.18,12.21,12.24,12.26,12.29,12.32,12.35,12.37,12.39,12.41,12.43,12.44,12.45,12.46,12.47,12.48,12.48,12.49,12.49,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.51,12.51,12.52,12.52,12.53,12.53,12.54,12.54,12.55,12.55,12.55,12.55,12.55,12.55,12.55,12.55,12.55,12.55,12.55,12.55,12.55,12.56,12.56,12.57,12.57,12.58,12.58,12.59,12.59,12.59,12.59,12.59,12.59,12.59,12.60,12.60,12.61,12.61,12.62,12.62,12.63,12.63,12.64,12.64,12.64,12.64,12.65,12.65,12.65,12.66,12.66,12.67,12.67,12.68,12.68,12.69,12.69,12.69,12.70,12.70,12.71,12.71,12.72,12.72,12.73,12.73,12.74,12.74,12.75,12.75,12.75,12.76,12.76,12.77,12.77,12.78,12.79,12.79,12.80,12.80,12.81,12.81,12.82,12.82,12.83,12.84,12.85,12.85,12.85,12.86,12.86,12.87,12.88,12.89,12.89,12.90,12.90,12.91,12.91,12.92,12.93,12.94,12.94,12.95,12.95,12.95,12.96,12.97,12.98,12.99,12.99,13.00,13.00,13.01,13.02,13.03,13.04,13.04,13.05,13.05,13.05,13.06,13.07,13.08,13.09,13.09,13.10,13.10,13.11,13.12,13.13,13.14,13.14,13.15,13.15,13.16,13.17,13.18,13.19,13.19,13.20,13.21,13.22,13.23,13.24,13.24,13.25,13.26,13.27,13.28,13.28,13.29,13.30,13.31,13.32,13.33,13.33,13.34,13.35,13.36,13.37,13.38,13.38,13.39,13.40,13.41,13.42,13.43,13.44,13.44,13.45,13.46,13.47,13.48,13.49,13.50,13.51,13.52,13.53,13.54,13.55,13.56,13.57,13.58,13.59,13.60,13.61,13.62,13.63,13.64,13.64,13.65,13.66,13.67,13.68,13.69,13.70,13.72,13.73,13.74,13.74,13.75,13.77,13.78,13.79,13.80,13.81,13.82,13.83,13.84,13.84,13.86,13.87,13.88,13.89,13.90,13.92,13.93,13.94,13.94,13.96,13.97,13.98,13.99,14.00,14.01,14.02,14.03,14.04,14.06,14.07,14.08,14.09,14.11,14.12,14.13,14.14,14.15,14.16,14.17,14.19,14.20,14.21,14.22,14.24,14.25,14.26,14.27,14.28,14.30,14.31,14.32,14.34,14.35,14.36,14.37,14.39,14.40,14.41,14.43,14.44,14.45,14.46,14.47,14.49,14.50,14.51,14.53,14.54,14.55,14.56,14.58,14.59,14.61,14.62,14.64,14.64,14.66,14.67,14.69,14.70,14.72,14.73,14.74,14.76,14.77,14.79,14.80,14.82,14.83,14.84,14.86,14.87,14.89,14.90,14.92,14.94,14.95,14.96,14.98,14.99,15.01,15.02,15.04,15.05,15.06,15.08,15.09,15.11,15.13,15.14,15.15,15.17,15.19,15.20,15.22,15.23,15.25,15.27,15.28,15.30,15.32,15.33,15.34,15.36,15.38,15.39,15.41,15.43,15.44,15.46,15.48,15.50,15.52,15.53,15.55,15.57,15.58,15.60,15.62,15.63,15.65,15.66,15.68,15.70,15.72,15.73,15.75,15.77,15.79,15.81,15.83,15.84,15.86,15.88,15.90,15.91,15.93,15.95,15.97,15.99,16.00,16.02,16.04,16.06,16.08,16.10,16.12,16.13,16.15,16.17,16.19,16.21,16.23,16.25,16.27,16.29,16.31,16.33,16.35,16.37,16.39,16.41,16.43,16.45,16.47,16.49,16.51,16.53,16.55,16.57,16.59,16.61,16.63,16.65,16.68,16.70,16.72,16.74,16.76,16.78,16.80,16.83,16.84,16.87,16.89,16.91,16.93,16.95,16.97,17.00,17.02,17.04,17.07,17.09,17.12,17.13,17.16,17.18,17.20,17.22,17.24,17.26,17.29,17.31,17.33,17.35,17.38,17.40,17.42,17.45,17.47,17.50,17.52,17.54,17.57,17.59,17.62,17.64,17.67,17.70,17.72,17.75,17.77,17.80,17.82,17.85,17.88,17.90,17.92,17.95,17.97,18.00,18.02,18.05,18.08,18.10,18.13,18.15,18.18,18.21,18.23,18.26,18.29,18.32,18.34,18.37,18.40,18.42,18.45,18.48,18.51,18.53,18.56,18.59,18.62,18.65,18.67,18.70,18.73,18.76,18.79,18.82,18.85,18.88,18.91,18.93,18.97,19.00,19.02,19.06,19.09,19.11,19.14,19.17,19.21,19.23,19.26,19.30,19.32,19.35,19.39,19.41,19.45,19.48,19.51,19.54,19.58,19.61,19.64,19.67,19.71,19.74,19.77,19.81,19.84,19.87,19.91,19.94,19.97,20.01,20.04,20.07,20.11,20.14,20.18,20.21,20.24,20.28,20.32,20.35,20.39,20.42,20.46,20.50,20.53,20.57,20.61,20.64,20.67,20.71,20.75,20.78,20.81,20.85,20.89,20.93,20.96,21.00,21.04,21.08,21.11,21.16,21.20,21.23,21.27,21.31,21.35,21.40,21.43,21.47,21.51,21.55,21.59,21.63,21.68,21.71,21.75,21.79,21.83,21.88,21.92,21.96,22.00,22.05,22.09,22.13,22.18,22.22,22.26,22.30,22.35,22.39,22.43,22.48,22.52,22.57,22.61,22.66,22.71,22.76,22.80,22.85,22.90,22.95,22.99,23.04,23.09,23.14,23.19,23.23,23.28,23.33,23.38,23.42,23.47,23.52,23.57,23.61,23.66,23.71,23.76,23.80,23.85,23.90,23.96,24.01,24.06,24.11,24.17,24.22,24.27,24.32,24.38,24.43,24.49,24.54,24.59,24.64,24.70,24.76,24.81,24.87,24.92,24.98,25.04,25.09,25.15,25.21,25.27,25.32,25.38,25.44,25.50,25.56,25.61,25.68,25.74,25.80,25.86,25.92,25.98,26.04,26.10,26.17,26.23,26.29,26.35,26.41,26.48,26.54,26.61,26.68,26.74,26.80,26.87,26.94,27.00,27.07,27.14,27.21,27.27,27.34,27.41,27.48,27.55,27.62,27.69,27.76,27.83,27.90,27.97,28.05,28.12,28.19,28.27,28.34,28.42,28.49,28.57,28.65,28.72,28.80,28.87,28.95,29.03,29.11,29.18,29.26,29.34,29.42,29.50,29.58,29.66,29.75,29.83,29.91,30.00,30.08,30.16,30.25,30.34,30.42,30.51,30.60,30.69,30.77,30.86,30.95,31.04,31.13,31.22,31.31,31.40,31.49,31.59,31.68,31.78,31.87,31.97,32.06,32.16,32.26,32.36,32.46,32.56,32.66,32.76,32.86,32.97,33.07,33.18,33.28,33.39,33.49,33.60,33.71,33.82,33.92,34.03,34.14,34.25,34.36,34.47,34.59,34.70,34.82,34.93,35.05,35.16,35.28,35.40,35.52,35.64,35.76,35.88,36.01,36.13,36.25,36.38,36.51,36.63,36.76,36.89,37.02,37.16,37.29,37.42,37.55,37.69,37.82,37.96,38.10,38.24,38.38,38.52,38.66,38.80,38.95,39.10,39.24,39.39,39.54,39.69,39.84,40.00,40.15,40.31,40.47,40.62,40.78,40.94,41.10,41.27,41.43,41.60,41.77,41.94,42.10,42.28,42.45,42.63,42.80,42.98,43.16,43.34,43.52,43.70,43.89,44.07,44.26,44.45,44.64,44.84,45.03,45.23,45.43,45.63,45.84,46.04,46.25,46.46,46.67,46.88,47.09,47.31,47.53,47.75,47.97,48.19,48.42,48.65,48.87,49.11,49.35,49.58,49.82,50.06,50.31,50.56,50.81,51.06,51.32,51.58,51.84,52.10,52.36,52.63,52.90,53.18,53.45,53.74,54.02,54.31,54.60,54.89,55.19,55.49,55.79,56.10,56.41,56.72,57.03,57.35,57.68,58.00,58.33,58.67,59.01,59.35,59.70,60.05,60.40,60.77,61.13,61.50,61.88,62.26,62.64,63.03,63.42,63.81,64.21,64.62,65.03,65.45,65.87,66.30,66.73,67.17,67.62,68.07,68.53,68.99,69.46,69.94,70.43,70.92,71.42,71.92,72.44,72.95,73.48,74.02,74.56,75.11,75.67,76.23,76.81,77.40,77.99,78.60,79.21,79.83,80.47,81.11,81.77,82.44,83.11,83.80,84.50,85.22,85.94,86.68,87.43,88.19,88.96,89.75,90.56,91.37,92.20,93.05,93.91,94.79,95.69,96.52,96.95,97.39,97.84,98.29];
                    
%     sec1 =  300;    sec2 = 800;     sec3 = 2048; 
%     M11 = 108;      M12 = 64;       M13 = 3; 
%                     M22 = 204;      M23 = 6; 
%                                     M33 = 9; 
%                 
%     sec1 =  593;    sec2 = 1075;    sec3 = 2048; 
%   	M11 = 54;       M12 = 32;       M13 = 3; 
%                     M22 = 102;      M23 = 6; 
%                                     M33 = 9; 
  
    if(0)   % Linear modulation scheme
%         sec1 =  593; sec2 = 1075; sec3 = 2048; 
%         M11 = 54;   M12 = 32;   M13 = 3; 
%                     M22 = 102;  M23 = 6; 
%                                 M33 = 9; 

        M3_list2 = round( [linspace(250,250,sec2) ...
                           linspace(250,M33,sec3-sec2)]);
        M2_list2 = round( [linspace(250,250,sec1) ...
                           linspace(250,M22,sec2-sec1) ...
                           linspace(M22,M23,sec3-sec2)]); 
        M1_list2 = round( [linspace(250,M11,sec1) ...
                           linspace(M11,M12,sec2-sec1) ...
                           linspace(M12,M13,sec3-sec2)]);

        disp(['M1 Length: ' num2str(length(M1_list2))]); 
        disp(['M2 Length: ' num2str(length(M2_list2))]); 
        disp(['M3 Length: ' num2str(length(M3_list2))]); 

        fprintf('slopes:  %6.4f   %6.4f   %6.4f\n', (250-M11)/(sec1), (M11-M12)/(sec2-sec1), (M12-M13)/(sec3-sec2)); 
        fprintf('slopes:  %6.4f   %6.4f   %6.4f\n', (250-250)/(sec1), (250-M22)/(sec2-sec1), (M22-M23)/(sec3-sec2)); 
        fprintf('slopes:  %6.4f   %6.4f   %6.4f\n', (250-250)/(sec1), (250-250)/(sec2-sec1), (250-M33)/(sec3-sec2)); 
        
        figure(11); hold off;
        plot(M1_list, 'lineWidth', 2); hold on;
        plot(M2_list, 'lineWidth', 2);
        plot(M3_list, 'lineWidth', 2);
        plot(M1_list2, 'k', 'lineWidth', 1);
        plot(M2_list2, 'lineWidth', 1);
        plot(M3_list2, 'lineWidth', 1);
        legend('M1', 'M2', 'M3', '', '', '');
        figure(11); hold off;
        plot(M1_list, 'lineWidth', 2); hold on;
        plot(M2_list, 'lineWidth', 2);
        plot(M3_list, 'lineWidth', 2);
        plot(M1_list2, 'k', 'lineWidth', 1);
        plot(M2_list2, 'k', 'lineWidth', 1);
        plot(M3_list2, 'k', 'lineWidth', 1);
        legend('M1', 'M2', 'M3');
    end
    
    VrecMin = 2.5; 
    VrecMax = 98.28; 
    
    ind = 2048 - bit; 
    
    M1 = M1_list(ind);
    M2 = M2_list(ind); 
    M3 = M3_list(ind);  
    
    nVs = (VrecFund(ind)-VrecMin)/VrecMax; 
    duty = ind/2048 * 100;
    
%         % THIS PLOTS THE FIGURE USED IN COMPEL 2020 PAPER
%         M1_list = [250,249,249,249,248,248,248,247,247,247,246,246,246,246,245,245,245,244,244,244,243,243,243,242,242,242,241,241,241,240,240,240,239,239,239,239,238,238,238,237,237,237,236,236,236,235,235,235,234,234,234,233,233,233,232,232,232,232,231,231,231,230,230,230,229,229,229,228,228,228,227,227,227,226,226,226,225,225,225,224,224,224,224,223,223,223,222,222,222,221,221,221,220,220,220,219,219,219,218,218,218,217,217,217,217,216,216,216,215,215,215,214,214,214,213,213,213,212,212,212,211,211,211,210,210,210,209,209,209,209,208,208,208,207,207,207,206,206,206,205,205,205,204,204,204,203,203,203,202,202,202,202,201,201,201,200,200,200,199,199,199,198,198,198,197,197,197,196,196,196,195,195,195,195,194,194,194,193,193,193,192,192,192,191,191,191,190,190,190,189,189,189,188,188,188,187,187,187,187,186,186,186,185,185,185,184,184,184,183,183,183,182,182,182,181,181,181,180,180,180,180,179,179,179,178,178,178,177,177,177,176,176,176,175,175,175,174,174,174,173,173,173,173,172,172,172,171,171,171,170,170,170,169,169,169,168,168,168,167,167,167,166,166,166,165,165,165,165,164,164,164,163,163,163,162,162,162,161,161,161,160,160,160,159,159,159,158,158,158,158,157,157,157,156,156,156,155,155,155,154,154,154,153,153,153,152,152,152,151,151,151,151,150,150,150,149,149,149,148,148,148,147,147,147,146,146,146,145,145,145,144,144,144,143,143,143,143,142,142,142,141,141,141,140,140,140,139,139,139,138,138,138,137,137,137,136,136,136,136,135,135,135,134,134,134,133,133,133,132,132,132,131,131,131,130,130,130,129,129,129,128,128,128,128,127,127,127,126,126,126,125,125,125,124,124,124,123,123,123,122,122,122,122,121,121,121,120,120,119,119,119,118,118,118,117,117,116,116,116,115,115,115,114,114,114,113,113,113,112,112,112,111,111,110,110,110,109,109,109,108,108,107,107,107,106,106,106,105,105,105,104,104,104,103,103,102,102,102,101,101,101,100,100,100,99,99,99,98,98,97,97,97,96,96,96,95,95,95,94,94,94,93,93,92,92,92,91,91,91,90,90,90,89,89,89,88,88,87,87,87,86,86,86,85,85,85,84,84,84,83,83,82,82,82,81,81,81,80,80,80,79,79,78,78,78,77,77,77,76,76,76,75,75,75,74,74,73,73,73,72,72,72,71,71,71,70,70,70,69,69,68,68,68,67,67,67,66,66,66,65,65,65,64,64,64,63,63,62,62,62,61,61,61,60,60,60,59,59,58,58,58,57,57,57,56,56,56,55,55,55,55,55,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,48,48,48,48,48,48,48,48,48,48,48,48,48,48,48,48,48,48,48,48,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,44,44,44,44,44,44,44,44,44,44,44,44,44,44,44,44,44,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,41,41,41,41,41,41,41,41,41,41,41,41,41,41,41,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,38,38,38,38,38,38,38,38,38,38,38,38,38,38,38,37,37,37,37,37,37,37,37,37,37,37,37,37,37,36,36,36,36,36,36,36,36,36,36,36,36,36,35,35,35,35,35,35,35,35,35,35,35,35,35,35,34,34,34,34,34,34,34,34,34,34,34,34,34,34,33,33,33,33,33,33,33,33,33,33,33,33,33,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,4,4,4,4,4,4,4,4,4,4,4,4,4,4];
%         M2_list = [250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,249,249,249,248,248,248,247,247,247,247,246,246,246,245,245,245,244,244,244,243,243,243,242,242,242,241,241,241,241,240,240,240,239,239,239,238,238,238,237,237,237,236,236,236,236,235,235,235,234,234,234,233,233,233,232,232,232,231,231,231,230,230,230,230,229,229,229,228,228,228,227,227,227,226,226,226,225,225,225,224,224,224,224,223,223,223,222,222,222,221,221,221,220,220,220,219,219,219,218,218,218,218,217,217,217,216,216,216,215,215,215,214,214,214,213,213,213,213,212,212,212,211,211,211,210,210,210,209,209,209,208,208,208,207,207,207,207,206,206,206,205,205,205,204,204,204,203,203,203,202,202,202,201,201,201,201,200,200,200,199,199,199,198,198,198,197,197,197,196,196,196,195,195,195,195,194,194,194,193,193,193,192,192,192,191,191,191,190,190,190,190,189,189,189,188,188,188,187,187,187,187,187,186,186,186,186,185,185,184,184,184,184,183,183,183,182,182,182,182,181,181,180,180,180,180,179,179,178,178,178,178,177,177,177,176,176,176,176,175,175,175,174,174,174,173,173,173,172,172,172,171,171,171,171,170,170,170,169,169,169,168,168,168,167,167,167,167,166,166,166,165,165,165,164,164,164,164,163,163,163,162,162,162,161,161,161,161,160,160,160,159,159,159,159,158,158,157,157,157,156,156,156,156,155,155,155,155,154,154,154,153,153,153,153,152,152,151,151,151,151,150,150,150,149,149,149,149,148,148,148,147,147,147,147,146,146,146,145,145,145,144,144,144,144,143,143,143,142,142,142,141,141,141,141,140,140,140,139,139,139,139,138,138,138,137,137,137,136,136,136,136,135,135,135,134,134,134,134,133,133,133,132,132,132,132,131,131,131,130,130,130,129,129,129,129,128,128,128,127,127,127,127,126,126,126,125,125,125,125,124,124,124,124,123,123,123,122,122,122,122,121,121,121,120,120,120,119,119,119,119,118,118,118,117,117,117,117,116,116,116,115,115,115,115,114,114,114,113,113,113,113,112,112,112,111,111,111,111,110,110,110,110,109,109,109,108,108,108,108,107,107,107,106,106,106,106,105,105,105,104,104,104,104,103,103,103,103,103,103,103,103,103,103,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,102,101,101,101,101,101,101,101,101,101,101,101,101,101,101,101,101,101,101,101,101,101,101,101,101,101,101,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,98,98,98,98,98,98,98,98,98,98,98,98,98,98,98,98,98,97,97,97,97,97,97,97,97,97,97,97,97,97,97,97,97,97,96,96,96,96,96,96,96,96,96,96,96,96,96,96,96,95,95,95,95,95,95,95,95,95,95,95,95,95,95,95,94,94,94,94,94,94,94,94,94,94,94,94,94,94,93,93,93,93,93,93,93,93,93,93,93,93,93,92,92,92,92,92,92,92,92,92,92,92,92,92,91,91,91,91,91,91,91,91,91,91,91,91,91,90,90,90,90,90,90,90,90,90,90,90,90,89,89,89,89,89,89,89,89,89,89,89,89,88,88,88,88,88,88,88,88,88,88,88,88,87,87,87,87,87,87,87,87,87,87,87,87,86,86,86,86,86,86,86,86,86,86,86,85,85,85,85,85,85,85,85,85,85,85,84,84,84,84,84,84,84,84,84,84,84,83,83,83,83,83,83,83,83,83,83,83,82,82,82,82,82,82,82,82,82,82,82,81,81,81,81,81,81,81,81,81,81,81,80,80,80,80,80,80,80,80,80,80,79,79,79,79,79,79,79,79,79,79,78,78,78,78,78,78,78,78,78,78,78,77,77,77,77,77,77,77,77,77,77,76,76,76,76,76,76,76,76,76,76,75,75,75,75,75,75,75,75,75,75,74,74,74,74,74,74,74,74,74,74,73,73,73,73,73,73,73,73,73,73,72,72,72,72,72,72,72,72,72,72,71,71,71,71,71,71,71,71,71,71,70,70,70,70,70,70,70,70,70,69,69,69,69,69,69,69,69,69,69,68,68,68,68,68,68,68,68,68,67,67,67,67,67,67,67,67,67,67,66,66,66,66,66,66,66,66,66,65,65,65,65,65,65,65,65,65,64,64,64,64,64,64,64,64,64,64,63,63,63,63,63,63,63,63,63,62,62,62,62,62,62,62,62,62,61,61,61,61,61,61,61,61,61,61,60,60,60,60,60,60,60,60,60,59,59,59,59,59,59,59,59,59,58,58,58,58,58,58,58,58,58,57,57,57,57,57,57,57,57,57,56,56,56,56,56,56,56,56,56,55,55,55,55,55,55,55,55,55,54,54,54,54,54,54,54,54,54,53,53,53,53,53,53,53,53,53,52,52,52,52,52,52,52,52,52,51,51,51,51,51,51,51,51,51,50,50,50,50,50,50,50,50,50,49,49,49,49,49,49,49,49,49,48,48,48,48,48,48,48,48,48,47,47,47,47,47,47,47,47,46,46,46,46,46,46,46,46,46,45,45,45,45,45,45,45,45,45,44,44,44,44,44,44,44,44,44,43,43,43,43,43,43,43,43,43,42,42,42,42,42,42,42,42,41,41,41,41,41,41,41,41,41,40,40,40,40,40,40,40,40,40,39,39,39,39,39,39,39,39,38,38,38,38,38,38,38,38,38,37,37,37,37,37,37,37,37,37,36,36,36,36,36,36,36,36,35,35,35,35,35,35,35,35,35,34,34,34,34,34,34,34,34,34,33,33,33,33,33,33,33,33,32,32,32,32,32,32,32,32,32,31,31,31,31,31,31,31,31,30,30,30,30,30,30,30,30,30,29,29,29,29,29,29,29,29,29,28,28,28,28,28,28,28,28,27,27,27,27,27,27,27,27,27,26,26,26,26,26,26,26,26,25,25,25,25,25,25,25,25,25,24,24,24,24,24,24,24,24,23,23,23,23,23,23,23,23,23,22,22,22,22,22,22,22,22,21,21,21,21,21,21,21,21,21,20,20,20,20,20,20,20,20,19,19,19,19,19,19,19,19,19,18,18,18,18,18,18,18,18,17,17,17,17,17,17,17,17,17,16,16,16,16,16,16,16,16,15,15,15,15,15,15,15,15,14,14,14,14,14,14,14,14,14,13,13,13,13,13,13,13,13,12,12,12,12,12,12];
%         M3_list = [250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,249,249,249,249,248,248,248,248,247,247,247,246,246,246,246,245,245,245,245,244,244,244,243,243,243,243,242,242,242,241,241,241,241,240,240,240,240,239,239,239,238,238,238,238,237,237,237,237,236,236,236,235,235,235,235,234,234,234,234,233,233,233,232,232,232,232,231,231,231,230,230,230,230,229,229,229,229,228,228,228,227,227,227,227,226,226,226,226,225,225,225,224,224,224,224,223,223,223,223,222,222,222,221,221,221,221,220,220,220,219,219,219,219,218,218,218,218,217,217,217,216,216,216,216,215,215,215,215,214,214,214,213,213,213,213,212,212,212,212,211,211,211,210,210,210,210,209,209,209,208,208,208,208,207,207,207,207,206,206,206,205,205,205,205,204,204,204,204,203,203,203,202,202,202,202,201,201,201,201,200,200,200,199,199,199,199,198,198,198,197,197,197,197,196,196,196,196,195,195,195,195,195,195,194,194,194,194,193,193,193,192,192,192,192,191,191,191,191,190,190,190,189,189,189,189,188,188,188,187,187,187,187,186,186,186,186,185,185,185,185,184,184,184,184,183,183,183,182,182,182,182,182,181,181,181,180,180,180,180,179,179,179,179,178,178,178,178,177,177,177,177,176,176,176,176,175,175,175,175,174,174,174,174,173,173,173,173,172,172,172,172,171,171,171,171,170,170,170,170,169,169,169,169,168,168,168,168,167,167,167,167,166,166,166,166,165,165,165,165,165,164,164,164,163,163,163,163,163,162,162,162,162,161,161,161,161,160,160,160,160,159,159,159,159,158,158,158,158,157,157,157,157,156,156,156,156,156,155,155,155,155,154,154,154,154,153,153,153,153,152,152,152,152,152,151,151,151,151,150,150,150,150,149,149,149,149,148,148,148,148,148,147,147,147,147,146,146,146,146,145,145,145,145,144,144,144,144,144,143,143,143,143,142,142,142,142,142,141,141,141,141,140,140,140,140,139,139,139,139,139,138,138,138,138,137,137,137,137,136,136,136,136,136,135,135,135,135,134,134,134,134,134,133,133,133,133,132,132,132,132,132,131,131,131,131,130,130,130,130,129,129,129,129,129,128,128,128,128,128,127,127,127,127,126,126,126,126,126,125,125,125,125,124,124,124,124,123,123,123,123,123,122,122,122,122,122,121,121,121,121,120,120,120,120,120,119,119,119,119,118,118,118,118,118,117,117,117,117,117,116,116,116,116,115,115,115,115,115,114,114,114,114,113,113,113,113,113,112,112,112,112,112,111,111,111,111,110,110,110,110,110,109,109,109,109,108,108,108,108,108,107,107,107,107,107,106,106,106,106,105,105,105,105,105,104,104,104,104,104,103,103,103,103,103,102,102,102,102,101,101,101,101,101,100,100,100,100,100,99,99,99,99,98,98,98,98,98,97,97,97,97,97,96,96,96,96,96,95,95,95,95,94,94,94,94,94,93,93,93,93,93,92,92,92,92,91,91,91,91,91,90,90,90,90,90,89,89,89,89,89,88,88,88,88,88,87,87,87,87,87,86,86,86,86,85,85,85,85,85,84,84,84,84,84,83,83,83,83,83,82,82,82,82,82,81,81,81,81,81,80,80,80,80,79,79,79,79,79,78,78,78,78,78,77,77,77,77,77,76,76,76,76,76,75,75,75,75,75,74,74,74,74,74,73,73,73,73,72,72,72,72,72,71,71,71,71,71,70,70,70,70,70,69,69,69,69,69,68,68,68,68,68,67,67,67,67,67,66,66,66,66,66,65,65,65,65,65,64,64,64,64,64,63,63,63,63,63,62,62,62,62,61,61,61,61,61,60,60,60,60,60,59,59,59,59,59,58,58,58,58,58,57,57,57,57,57,56,56,56,56,56,55,55,55,55,55,54,54,54,54,54,53,53,53,53,53,52,52,52,52,52,51,51,51,51,51,50,50,50,50,50,49,49,49,49,49,48,48,48,48,48,47,47,47,47,47,46,46,46,46,46,45,45,45,45,45,44,44,44,44,44,43,43,43,43,43,42,42,42,42,42,41,41,41,41,41,40,40,40,40,40,39,39,39,39,39,38,38,38,38,38,37,37,37,37,37,36,36,36,36,36,35,35,35,35,35,34,34,34,34,34,33,33,33,33,33,32,32,32,32,32,31,31,31,31,31,30,30,30,30,30,29,29,29,29,29,28,28,28,28,28,27,27,27,27,27,26,26,26,26,26,25,25,25,25,25,24,24,24,24,24,23,23,23,23,23,22,22,22,22,22,21,21,21,21,21,21,20,20,20];VrecFund = [2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.51, 2.51, 2.52, 2.52, 2.53, 2.53, 2.54, 2.54, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.55, 2.56, 2.56, 2.57, 2.57, 2.58, 2.58, 2.59, 2.59, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.60, 2.61, 2.61, 2.62, 2.62, 2.63, 2.63, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.64, 2.65, 2.65, 2.66, 2.66, 2.67, 2.67, 2.68, 2.68, 2.69, 2.69, 2.69, 2.69, 2.69, 2.69, 2.69, 2.69, 2.69, 2.69, 2.69, 2.70, 2.70, 2.70, 2.71, 2.71, 2.72, 2.72, 2.73, 2.73, 2.74, 2.74, 2.74, 2.74, 2.74, 2.74, 2.74, 2.74, 2.74, 2.75, 2.75, 2.76, 2.76, 2.77, 2.77, 2.78, 2.78, 2.79, 2.79, 2.79, 2.79, 2.79, 2.79, 2.79, 2.80, 2.80, 2.80, 2.81, 2.81, 2.82, 2.82, 2.83, 2.83, 2.83, 2.83, 2.83, 2.84, 2.84, 2.85, 2.85, 2.86, 2.86, 2.87, 2.87, 2.88, 2.88, 2.88, 2.88, 2.88, 2.89, 2.89, 2.90, 2.90, 2.90, 2.91, 2.91, 2.92, 2.92, 2.93, 2.93, 2.93, 2.93, 2.94, 2.94, 2.95, 2.95, 2.96, 2.96, 2.97, 2.97, 2.98, 2.98, 2.98, 2.99, 2.99, 3.00, 3.00, 3.00, 3.01, 3.01, 3.02, 3.02, 3.03, 3.03, 3.04, 3.04, 3.05, 3.05, 3.06, 3.06, 3.07, 3.08, 3.08, 3.09, 3.09, 3.10, 3.10, 3.10, 3.11, 3.11, 3.12, 3.13, 3.13, 3.14, 3.14, 3.15, 3.15, 3.16, 3.16, 3.17, 3.18, 3.18, 3.19, 3.19, 3.20, 3.20, 3.20, 3.21, 3.22, 3.23, 3.23, 3.24, 3.24, 3.25, 3.26, 3.27, 3.28, 3.28, 3.29, 3.29, 3.30, 3.30, 3.31, 3.32, 3.32, 3.33, 3.33, 3.34, 3.34, 3.35, 3.36, 3.37, 3.38, 3.38, 3.39, 3.40, 3.40, 3.41, 3.42, 3.43, 3.43, 3.44, 3.45, 3.46, 3.47, 3.47, 3.48, 3.48, 3.49, 3.50, 3.51, 3.52, 3.52, 3.53, 3.54, 3.55, 3.56, 3.57, 3.57, 3.58, 3.59, 3.59, 3.60, 3.61, 3.62, 3.63, 3.64, 3.65, 3.66, 3.67, 3.68, 3.69, 3.69, 3.70, 3.71, 3.72, 3.73, 3.74, 3.75, 3.76, 3.77, 3.78, 3.79, 3.80, 3.81, 3.82, 3.83, 3.84, 3.85, 3.86, 3.87, 3.88, 3.89, 3.89, 3.91, 3.92, 3.93, 3.94, 3.95, 3.97, 3.98, 3.99, 4.00, 4.01, 4.02, 4.03, 4.05, 4.06, 4.07, 4.08, 4.09, 4.11, 4.12, 4.13, 4.14, 4.16, 4.17, 4.18, 4.19, 4.20, 4.21, 4.23, 4.24, 4.26, 4.27, 4.28, 4.29, 4.31, 4.32, 4.34, 4.35, 4.36, 4.38, 4.39, 4.40, 4.42, 4.43, 4.45, 4.46, 4.48, 4.49, 4.51, 4.52, 4.54, 4.56, 4.57, 4.59, 4.60, 4.62, 4.63, 4.65, 4.67, 4.68, 4.70, 4.71, 4.73, 4.75, 4.77, 4.79, 4.80, 4.82, 4.84, 4.86, 4.88, 4.89, 4.91, 4.93, 4.95, 4.97, 4.99, 5.01, 5.03, 5.05, 5.07, 5.09, 5.11, 5.13, 5.15, 5.17, 5.19, 5.21, 5.24, 5.26, 5.28, 5.30, 5.32, 5.35, 5.37, 5.39, 5.41, 5.43, 5.46, 5.48, 5.51, 5.53, 5.56, 5.58, 5.61, 5.64, 5.67, 5.69, 5.71, 5.74, 5.77, 5.79, 5.82, 5.85, 5.88, 5.90, 5.93, 5.97, 5.99, 6.02, 6.05, 6.08, 6.11, 6.14, 6.17, 6.20, 6.24, 6.27, 6.30, 6.33, 6.37, 6.40, 6.43, 6.47, 6.51, 6.54, 6.58, 6.61, 6.65, 6.69, 6.73, 6.77, 6.80, 6.84, 6.88, 6.92, 6.97, 7.00, 7.04, 7.08, 7.12, 7.17, 7.21, 7.25, 7.29, 7.32, 7.35, 7.38, 7.40, 7.43, 7.45, 7.46, 7.47, 7.47, 7.48, 7.48, 7.49, 7.49, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.51, 7.51, 7.52, 7.52, 7.53, 7.53, 7.54, 7.54, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.55, 7.56, 7.56, 7.57, 7.57, 7.57, 7.58, 7.58, 7.59, 7.59, 7.59, 7.59, 7.59, 7.59, 7.59, 7.59, 7.59, 7.60, 7.60, 7.61, 7.61, 7.62, 7.62, 7.63, 7.63, 7.64, 7.64, 7.64, 7.64, 7.64, 7.64, 7.64, 7.65, 7.65, 7.66, 7.66, 7.67, 7.67, 7.67, 7.68, 7.68, 7.69, 7.69, 7.69, 7.69, 7.70, 7.70, 7.71, 7.71, 7.72, 7.72, 7.73, 7.73, 7.74, 7.74, 7.74, 7.74, 7.75, 7.75, 7.76, 7.76, 7.77, 7.77, 7.77, 7.78, 7.78, 7.79, 7.79, 7.80, 7.80, 7.81, 7.81, 7.82, 7.82, 7.83, 7.83, 7.84, 7.84, 7.85, 7.85, 7.86, 7.86, 7.87, 7.87, 7.87, 7.88, 7.89, 7.89, 7.90, 7.90, 7.91, 7.91, 7.92, 7.92, 7.93, 7.94, 7.94, 7.95, 7.95, 7.96, 7.96, 7.97, 7.97, 7.98, 7.99, 7.99, 8.00, 8.00, 8.01, 8.01, 8.02, 8.03, 8.04, 8.04, 8.05, 8.05, 8.06, 8.07, 8.07, 8.08, 8.09, 8.09, 8.10, 8.10, 8.11, 8.12, 8.13, 8.14, 8.14, 8.15, 8.15, 8.16, 8.17, 8.18, 8.18, 8.19, 8.19, 8.20, 8.21, 8.22, 8.23, 8.23, 8.24, 8.24, 8.25, 8.26, 8.27, 8.28, 8.28, 8.29, 8.30, 8.31, 8.32, 8.33, 8.33, 8.34, 8.35, 8.36, 8.37, 8.37, 8.38, 8.38, 8.39, 8.40, 8.41, 8.42, 8.43, 8.44, 8.45, 8.46, 8.47, 8.47, 8.48, 8.49, 8.50, 8.51, 8.52, 8.52, 8.53, 8.54, 8.55, 8.56, 8.57, 8.58, 8.59, 8.60, 8.61, 8.62, 8.63, 8.64, 8.65, 8.66, 8.67, 8.67, 8.68, 8.69, 8.70, 8.71, 8.72, 8.73, 8.74, 8.75, 8.76, 8.77, 8.78, 8.79, 8.80, 8.81, 8.82, 8.83, 8.84, 8.85, 8.87, 8.87, 8.88, 8.89, 8.91, 8.92, 8.93, 8.94, 8.95, 8.97, 8.97, 8.98, 9.00, 9.01, 9.02, 9.03, 9.04, 9.06, 9.07, 9.07, 9.08, 9.10, 9.11, 9.12, 9.13, 9.14, 9.16, 9.17, 9.18, 9.19, 9.21, 9.22, 9.23, 9.25, 9.26, 9.27, 9.28, 9.30, 9.31, 9.32, 9.33, 9.35, 9.36, 9.37, 9.38, 9.40, 9.41, 9.42, 9.44, 9.45, 9.46, 9.48, 9.49, 9.50, 9.52, 9.53, 9.55, 9.56, 9.57, 9.59, 9.60, 9.62, 9.63, 9.65, 9.66, 9.67, 9.69, 9.70, 9.72, 9.73, 9.75, 9.76, 9.77, 9.79, 9.81, 9.82, 9.84, 9.86, 9.87, 9.89, 9.90, 9.92, 9.94, 9.95, 9.96, 9.98,10.00,10.01,10.03,10.05,10.06,10.08,10.09,10.11,10.13,10.14,10.16,10.18,10.19,10.21,10.23,10.25,10.26,10.28,10.30,10.32,10.34,10.36,10.37,10.39,10.41,10.43,10.45,10.46,10.48,10.50,10.52,10.54,10.56,10.57,10.59,10.61,10.63,10.65,10.67,10.69,10.71,10.73,10.75,10.77,10.79,10.81,10.83,10.85,10.87,10.89,10.91,10.93,10.96,10.97,10.99,11.02,11.04,11.06,11.08,11.10,11.12,11.15,11.16,11.19,11.21,11.24,11.26,11.28,11.31,11.33,11.35,11.38,11.40,11.43,11.45,11.47,11.50,11.52,11.55,11.57,11.59,11.62,11.64,11.66,11.69,11.71,11.74,11.76,11.78,11.81,11.83,11.85,11.88,11.91,11.94,11.96,11.99,12.01,12.04,12.07,12.10,12.13,12.15,12.18,12.21,12.24,12.26,12.29,12.32,12.35,12.37,12.39,12.41,12.43,12.44,12.45,12.46,12.47,12.48,12.48,12.49,12.49,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.50,12.51,12.51,12.52,12.52,12.53,12.53,12.54,12.54,12.55,12.55,12.55,12.55,12.55,12.55,12.55,12.55,12.55,12.55,12.55,12.55,12.55,12.56,12.56,12.57,12.57,12.58,12.58,12.59,12.59,12.59,12.59,12.59,12.59,12.59,12.60,12.60,12.61,12.61,12.62,12.62,12.63,12.63,12.64,12.64,12.64,12.64,12.65,12.65,12.65,12.66,12.66,12.67,12.67,12.68,12.68,12.69,12.69,12.69,12.70,12.70,12.71,12.71,12.72,12.72,12.73,12.73,12.74,12.74,12.75,12.75,12.75,12.76,12.76,12.77,12.77,12.78,12.79,12.79,12.80,12.80,12.81,12.81,12.82,12.82,12.83,12.84,12.85,12.85,12.85,12.86,12.86,12.87,12.88,12.89,12.89,12.90,12.90,12.91,12.91,12.92,12.93,12.94,12.94,12.95,12.95,12.95,12.96,12.97,12.98,12.99,12.99,13.00,13.00,13.01,13.02,13.03,13.04,13.04,13.05,13.05,13.05,13.06,13.07,13.08,13.09,13.09,13.10,13.10,13.11,13.12,13.13,13.14,13.14,13.15,13.15,13.16,13.17,13.18,13.19,13.19,13.20,13.21,13.22,13.23,13.24,13.24,13.25,13.26,13.27,13.28,13.28,13.29,13.30,13.31,13.32,13.33,13.33,13.34,13.35,13.36,13.37,13.38,13.38,13.39,13.40,13.41,13.42,13.43,13.44,13.44,13.45,13.46,13.47,13.48,13.49,13.50,13.51,13.52,13.53,13.54,13.55,13.56,13.57,13.58,13.59,13.60,13.61,13.62,13.63,13.64,13.64,13.65,13.66,13.67,13.68,13.69,13.70,13.72,13.73,13.74,13.74,13.75,13.77,13.78,13.79,13.80,13.81,13.82,13.83,13.84,13.84,13.86,13.87,13.88,13.89,13.90,13.92,13.93,13.94,13.94,13.96,13.97,13.98,13.99,14.00,14.01,14.02,14.03,14.04,14.06,14.07,14.08,14.09,14.11,14.12,14.13,14.14,14.15,14.16,14.17,14.19,14.20,14.21,14.22,14.24,14.25,14.26,14.27,14.28,14.30,14.31,14.32,14.34,14.35,14.36,14.37,14.39,14.40,14.41,14.43,14.44,14.45,14.46,14.47,14.49,14.50,14.51,14.53,14.54,14.55,14.56,14.58,14.59,14.61,14.62,14.64,14.64,14.66,14.67,14.69,14.70,14.72,14.73,14.74,14.76,14.77,14.79,14.80,14.82,14.83,14.84,14.86,14.87,14.89,14.90,14.92,14.94,14.95,14.96,14.98,14.99,15.01,15.02,15.04,15.05,15.06,15.08,15.09,15.11,15.13,15.14,15.15,15.17,15.19,15.20,15.22,15.23,15.25,15.27,15.28,15.30,15.32,15.33,15.34,15.36,15.38,15.39,15.41,15.43,15.44,15.46,15.48,15.50,15.52,15.53,15.55,15.57,15.58,15.60,15.62,15.63,15.65,15.66,15.68,15.70,15.72,15.73,15.75,15.77,15.79,15.81,15.83,15.84,15.86,15.88,15.90,15.91,15.93,15.95,15.97,15.99,16.00,16.02,16.04,16.06,16.08,16.10,16.12,16.13,16.15,16.17,16.19,16.21,16.23,16.25,16.27,16.29,16.31,16.33,16.35,16.37,16.39,16.41,16.43,16.45,16.47,16.49,16.51,16.53,16.55,16.57,16.59,16.61,16.63,16.65,16.68,16.70,16.72,16.74,16.76,16.78,16.80,16.83,16.84,16.87,16.89,16.91,16.93,16.95,16.97,17.00,17.02,17.04,17.07,17.09,17.12,17.13,17.16,17.18,17.20,17.22,17.24,17.26,17.29,17.31,17.33,17.35,17.38,17.40,17.42,17.45,17.47,17.50,17.52,17.54,17.57,17.59,17.62,17.64,17.67,17.70,17.72,17.75,17.77,17.80,17.82,17.85,17.88,17.90,17.92,17.95,17.97,18.00,18.02,18.05,18.08,18.10,18.13,18.15,18.18,18.21,18.23,18.26,18.29,18.32,18.34,18.37,18.40,18.42,18.45,18.48,18.51,18.53,18.56,18.59,18.62,18.65,18.67,18.70,18.73,18.76,18.79,18.82,18.85,18.88,18.91,18.93,18.97,19.00,19.02,19.06,19.09,19.11,19.14,19.17,19.21,19.23,19.26,19.30,19.32,19.35,19.39,19.41,19.45,19.48,19.51,19.54,19.58,19.61,19.64,19.67,19.71,19.74,19.77,19.81,19.84,19.87,19.91,19.94,19.97,20.01,20.04,20.07,20.11,20.14,20.18,20.21,20.24,20.28,20.32,20.35,20.39,20.42,20.46,20.50,20.53,20.57,20.61,20.64,20.67,20.71,20.75,20.78,20.81,20.85,20.89,20.93,20.96,21.00,21.04,21.08,21.11,21.16,21.20,21.23,21.27,21.31,21.35,21.40,21.43,21.47,21.51,21.55,21.59,21.63,21.68,21.71,21.75,21.79,21.83,21.88,21.92,21.96,22.00,22.05,22.09,22.13,22.18,22.22,22.26,22.30,22.35,22.39,22.43,22.48,22.52,22.57,22.61,22.66,22.71,22.76,22.80,22.85,22.90,22.95,22.99,23.04,23.09,23.14,23.19,23.23,23.28,23.33,23.38,23.42,23.47,23.52,23.57,23.61,23.66,23.71,23.76,23.80,23.85,23.90,23.96,24.01,24.06,24.11,24.17,24.22,24.27,24.32,24.38,24.43,24.49,24.54,24.59,24.64,24.70,24.76,24.81,24.87,24.92,24.98,25.04,25.09,25.15,25.21,25.27,25.32,25.38,25.44,25.50,25.56,25.61,25.68,25.74,25.80,25.86,25.92,25.98,26.04,26.10,26.17,26.23,26.29,26.35,26.41,26.48,26.54,26.61,26.68,26.74,26.80,26.87,26.94,27.00,27.07,27.14,27.21,27.27,27.34,27.41,27.48,27.55,27.62,27.69,27.76,27.83,27.90,27.97,28.05,28.12,28.19,28.27,28.34,28.42,28.49,28.57,28.65,28.72,28.80,28.87,28.95,29.03,29.11,29.18,29.26,29.34,29.42,29.50,29.58,29.66,29.75,29.83,29.91,30.00,30.08,30.16,30.25,30.34,30.42,30.51,30.60,30.69,30.77,30.86,30.95,31.04,31.13,31.22,31.31,31.40,31.49,31.59,31.68,31.78,31.87,31.97,32.06,32.16,32.26,32.36,32.46,32.56,32.66,32.76,32.86,32.97,33.07,33.18,33.28,33.39,33.49,33.60,33.71,33.82,33.92,34.03,34.14,34.25,34.36,34.47,34.59,34.70,34.82,34.93,35.05,35.16,35.28,35.40,35.52,35.64,35.76,35.88,36.01,36.13,36.25,36.38,36.51,36.63,36.76,36.89,37.02,37.16,37.29,37.42,37.55,37.69,37.82,37.96,38.10,38.24,38.38,38.52,38.66,38.80,38.95,39.10,39.24,39.39,39.54,39.69,39.84,40.00,40.15,40.31,40.47,40.62,40.78,40.94,41.10,41.27,41.43,41.60,41.77,41.94,42.10,42.28,42.45,42.63,42.80,42.98,43.16,43.34,43.52,43.70,43.89,44.07,44.26,44.45,44.64,44.84,45.03,45.23,45.43,45.63,45.84,46.04,46.25,46.46,46.67,46.88,47.09,47.31,47.53,47.75,47.97,48.19,48.42,48.65,48.87,49.11,49.35,49.58,49.82,50.06,50.31,50.56,50.81,51.06,51.32,51.58,51.84,52.10,52.36,52.63,52.90,53.18,53.45,53.74,54.02,54.31,54.60,54.89,55.19,55.49,55.79,56.10,56.41,56.72,57.03,57.35,57.68,58.00,58.33,58.67,59.01,59.35,59.70,60.05,60.40,60.77,61.13,61.50,61.88,62.26,62.64,63.03,63.42,63.81,64.21,64.62,65.03,65.45,65.87,66.30,66.73,67.17,67.62,68.07,68.53,68.99,69.46,69.94,70.43,70.92,71.42,71.92,72.44,72.95,73.48,74.02,74.56,75.11,75.67,76.23,76.81,77.40,77.99,78.60,79.21,79.83,80.47,81.11,81.77,82.44,83.11,83.80,84.50,85.22,85.94,86.68,87.43,88.19,88.96,89.75,90.56,91.37,92.20,93.05,93.91,94.79,95.69,96.52,96.95,97.39,97.84,98.29];
% 
%         M1plot = abs(M1_list-250)./250*Ts/2*1e6 % *100/Ts*2/1e6;
%         M2plot = abs(M2_list-250)./250*Ts/2*1e6 % *100/Ts*2/1e6;
%         M3plot = abs(M3_list-250)./250*Ts/2*1e6 % *100/Ts*2/1e6;
% 
%         KM3(1) = 0.3333;	KM2(1) = 0;         KM1(1) = 0;
%         KM3(2) = 0.0485;	KM2(2) = 0.3101;    KM1(2) = 0; 
%         KM3(3) = 0.0284; 	KM2(3) = 0.0923;  	KM1(3) = 0.2333; 
% 
%         high3 = [250,54,32,4];          high3 = abs(high3-250)./250*Ts/2*1e6 % *100/Ts*2/1e6; 
%         high2 = [250,250,103,12];       high2 = abs(high2-250)./250*Ts/2*1e6 % *100/Ts*2/1e6; 
%         high1 = [250,250,20];           high1 = abs(high1-250)./250*Ts/2*1e6 % *100/Ts*2/1e6; 
% 
%         breaks3 = [0,593,1063,2048];
%         breaks2 = [0,581,1065,2048]; 
%         breaks1 = [0,1048,2048];
% 
%         ff = figure(1); hold off; font = 'Times'; fntSize = 15;
%         set(ff, 'DefaultTextFontName', font, 'DefaultAxesFontName', font);
%         hh = plot(M1plot, 'lineWidth', 6, 'color', [210,210,210]./255); hold on; 
%         plot(M2plot, 'lineWidth', 6, 'color', [210,210,210]./255);
%         plot(M3plot, 'lineWidth', 6, 'color', [210,210,210]./255);
%         plot(breaks3, high3, '--k', 'lineWidth', 2); 
%             text(330, 27./100*Ts/2*1e6, 'K_{M1}', 'fontSize', fntSize-1);
%             x = 300;   y1 = 40;  y2 = 20; plot([x,x], [y1,y2]./100*Ts/2*1e6, 'k'); 
%             x1 = 150; x2 = 300;  y = 20; plot([x1,x2], [y,y]./100*Ts/2*1e6, 'k'); 
%         plot(breaks2, high2, '--k', 'lineWidth', 2); 
%             text(890, 22./100*Ts/2*1e6, 'K_{M2}', 'fontSize', fntSize-1);
%             x = 860;   y1 = 33;  y2 = 15; plot([x,x], [y1,y2]./100*Ts/2*1e6, 'k'); 
%             x1 = 710; x2 = 860;  y = 15; plot([x1,x2], [y,y]./100*Ts/2*1e6, 'k'); 
%         plot(breaks1, high1, '--k', 'lineWidth', 2); 
%             text(1350, 19./100*Ts/2*1e6, 'K_{M3}', 'fontSize', fntSize-1);
%             x = 1320;   y1 = 25;  y2 = 12; plot([x,x], [y1,y2]./100*Ts/2*1e6, 'k'); 
%             x1 = 1170; x2 = 1320;  y = 12; plot([x1,x2], [y,y]./100*Ts/2*1e6, 'k'); 
%         % text(480, 1.1*Ts/4*1e6, 'M1', 'fontSize', fntSize-1);
%         % text(960, 0.8*Ts/4*1e6, 'M2', 'fontSize', fntSize-1);
%         % text(1385, 0.52*Ts/4*1e6, 'M3', 'fontSize', fntSize-1);
%         text(360, 70./100*Ts/2*1e6, 'M1', 'fontSize', fntSize-1);
%         text(850, 55./100*Ts/2*1e6, 'M2', 'fontSize', fntSize-1);
%         text(1280, 43./100*Ts/2*1e6, 'M3', 'fontSize', fntSize-1);
%         xlim([-10, 2058]); 
%         ylim([-Ts/128*1e6, Ts/2*1e6]); 
%         xlabel('11 Bit Modulator Input'); 
%         ylabel('On-Time  [{\mu}s]');          
%         set(gca, 'Layer', 'Top');
%             ax = ancestor(hh, 'axes'); 
%             yticks([0, Ts/8*1e6, Ts/4*1e6, 3*Ts/8*1e6, Ts/2*1e6])
%             yticklabels({'0', '0.83', '1.67', '2.50', 'P/2'})
%         %     yticks([0, 25, 50, 75, 100])
%             xrule = ax.XAxis; xrule.FontSize = fntSize-1; 
%             yrule = ax.YAxis; yrule.FontSize = fntSize;  

    
end
% --------------------------------------------------------------


% % --------------- Functions (Vrec Referenced) ------------------
% % -------------------- 7 Level Intervals -----------------------
% function [td,te,te2,topA,topB,inPair0,inPolD,inPolP,Vrs,Us] = ...
%             modInts7(clk2,Pclk150,ph,M1,M2,M3,Az,Ap1,Ap2,Ap3,An1,An2,An3,B)
% 
%     if(ph==M1 || ph==M2 || ph==M3)
%         error(['Interval Definitions: edge times are invalid: ph=' num2str(ph)...
%             ', M1=' num2str(M1) ', M2=' num2str(M2) ', M3=' num2str(M3)]); 
%     elseif(abs(ph)<M1)
%         if(ph>0)
%             td = [ph,M1-ph,M2-M1,M3-M2,clk2-2*M3,M3-M2,M2-M1,M1]; td = [td td];
%             topA = cat(3,Az,Az,Ap1,Ap2,Ap3,Ap2,Ap1,Az,      Az,Az,An1,An2,An3,An2,An1,Az); 
%             topB = cat(3,-B,B,B,B,B,B,B,B,                  B,-B,-B,-B,-B,-B,-B,-B); 
%             Vrs = [0,0,1,2,3,2,1,0,                         0,0,-1,-2,-3,-2,-1,0];
%             Us = [-1,1,1,1,1,1,1,1,                         1,-1,-1,-1,-1,-1,-1,-1]; 
%             inPair0 = [2,3,4,5,6,7; 3,4,5,6,7,8]; 
%                 inPair02 = inPair0 + 8; 
%                 inPair0 = [inPair0 inPair02]';
%             inPolD = [1,1,1,0,0,0,                          1,1,1,0,0,0]'; 
%             inPolP = [0,0,0,0,0,0,                          0,0,0,0,0,0]';
%         else
%             td = [ph,M1-ph,M2-M1,M3-M2,clk2-2*M3,M3-M2,M2-M1,M1]; td = [td td];
%             topA = cat(3,Az,Az,Ap1,Ap2,Ap3,Ap2,Ap1,Az,      Az,Az,An1,An2,An3,An2,An1,Az); 
%             topB = cat(3,-B,B,B,B,B,B,B,B,                  B,-B,-B,-B,-B,-B,-B,-B); 
%             Vrs = [0,0,1,2,3,2,1,0,                         0,0,-1,-2,-3,-2,-1,0];
%             Us = [-1,1,1,1,1,1,1,1,                         1,-1,-1,-1,-1,-1,-1,-1]; 
%             inPair0 = [2,3,4,5,6,7; 3,4,5,6,7,8]; 
%                 inPair02 = inPair0 + 8; 
%                 inPair0 = [inPair0 inPair02]';
%             inPolD = [1,1,1,0,0,0,                          1,1,1,0,0,0]'; 
%             inPolP = [0,0,0,0,0,0,                          0,0,0,0,0,0]';
%         end
%     elseif(abs(ph)<M2) 
%         if(ph>0) 
%             td = [M1,ph-M1,M2-ph,M3-M2,clk2-2*M3,M3-M2,M2-M1,M1]; td = [td td];
%             topA = cat(3,Az,Ap1,Ap1,Ap2,Ap3,Ap2,Ap1,Az,     Az,An1,An1,An2,An3,An2,An1,Az); 
%             topB = cat(3,-B,-B,B,B,B,B,B,B,                 B,B,-B,-B,-B,-B,-B,-B); 
%             Vrs = [0,1,1,2,3,2,1,0,                         0,-1,-1,-2,-3,-2,-1,0];
%             Us = [-1,-1,1,1,1,1,1,1,                        1,1,-1,-1,-1,-1,-1,-1]; 
%             inPair0 = [1,3,4,5,6,7; 2,4,5,6,7,8];             
%                 inPair02 = inPair0 + 8; 
%                 inPair0 = [inPair0 inPair02]';
%             inPolD = [1,1,1,0,0,0,                          1,1,1,0,0,0]'; 
%             inPolP = [0,0,0,0,0,0,                          0,0,0,0,0,0]'; 
%         else
%             td = [M1,ph-M1,M2-ph,M3-M2,clk2-2*M3,M3-M2,M2-M1,M1]; td = [td td];
%             topA = cat(3,Az,Ap1,Ap1,Ap2,Ap3,Ap2,Ap1,Az,     Az,An1,An1,An2,An3,An2,An1,Az); 
%             topB = cat(3,-B,-B,B,B,B,B,B,B,                 B,B,-B,-B,-B,-B,-B,-B); 
%             Vrs = [0,1,1,2,3,2,1,0,                         0,-1,-1,-2,-3,-2,-1,0];
%             Us = [-1,-1,1,1,1,1,1,1,                        1,1,-1,-1,-1,-1,-1,-1]; 
%             inPair0 = [1,3,4,5,6,7; 2,4,5,6,7,8];             
%                 inPair02 = inPair0 + 8; 
%                 inPair0 = [inPair0 inPair02]';
%             inPolD = [1,1,1,0,0,0,                          1,1,1,0,0,0]'; 
%             inPolP = [0,0,0,0,0,0,                          0,0,0,0,0,0]';
%         end
%     elseif(abs(ph)<M3)
%         if(ph>0)
%             td = [M1,M2-M1,ph-M2,M3-ph,clk2-2*M3,M3-M2,M2-M1,M1]; td = [td td];
%             topA = cat(3,Az,Ap1,Ap2,Ap2,Ap3,Ap2,Ap1,Az,     Az,An1,An2,An2,An3,An2,An1,Az); 
%             topB = cat(3,-B,-B,-B,B,B,B,B,B,                B,B,B,-B,-B,-B,-B,-B); 
%             Vrs = [0,1,2,2,3,2,1,0,                         0,-1,-2,-2,-3,-2,-1,0];
%             Us = [-1,-1,-1,1,1,1,1,1,                       1,1,1,-1,-1,-1,-1,-1]; 
%             inPair0 = [1,2,4,5,6,7; 2,3,5,6,7,8];           
%                 inPair02 = inPair0 + 8; 
%                 inPair0 = [inPair0 inPair02]';
%             inPolD = [1,1,1,0,0,0,                          1,1,1,0,0,0]'; 
%             inPolP = [0,0,0,0,0,0,                          0,0,0,0,0,0]';
%         else
%             td = [M1,M2-M1,ph-M2,M3-ph,clk2-2*M3,M3-M2,M2-M1,M1]; td = [td td];
%             topA = cat(3,Az,Ap1,Ap2,Ap2,Ap3,Ap2,Ap1,Az,     Az,An1,An2,An2,An3,An2,An1,Az); 
%             topB = cat(3,-B,-B,-B,B,B,B,B,B,                B,B,B,-B,-B,-B,-B,-B); 
%             Vrs = [0,1,2,2,3,2,1,0,                         0,-1,-2,-2,-3,-2,-1,0];
%             Us = [-1,-1,-1,1,1,1,1,1,                       1,1,1,-1,-1,-1,-1,-1]; 
%             inPair0 = [1,2,4,5,6,7; 2,3,5,6,7,8];           
%                 inPair02 = inPair0 + 8; 
%                 inPair0 = [inPair0 inPair02]';
%             inPolD = [1,1,1,0,0,0,                          1,1,1,0,0,0]'; 
%             inPolP = [0,0,0,0,0,0,                          0,0,0,0,0,0]';
%         end
%     elseif(abs(ph)>M3)
%         if(ph>0)
%             td = [M1,M2-M1,M3-M2,ph-M3,clk2-M3-ph,M3-M2,M2-M1,M1]; td = [td td];
%             topA = cat(3,Az,Ap1,Ap2,Ap3,Ap3,Ap2,Ap1,Az,     Az,An1,An2,An3,An3,An2,An1,Az); 
%             topB = cat(3,-B,-B,-B,-B,B,B,B,B,               B,B,B,B,-B,-B,-B,-B); 
%             Vrs = [0,1,2,3,3,2,1,0,                         0,-1,-2,-3,-3,-2,-1,0];
%             Us = [-1,-1,-1,-1,1,1,1,1,                      1,1,1,1,-1,-1,-1,-1]; 
%             inPair0 = [1,2,3,5,6,7; 2,3,4,6,7,8];   
%                 inPair02 = inPair0 + 8; inPair0 = [inPair0 inPair02]';     
%             inPolD = [1,1,1,0,0,0,                          1,1,1,0,0,0]'; 
%             inPolP = [0,0,0,0,0,0,                          0,0,0,0,0,0]';
%         else
%             td = [M1,M2-M1,M3-M2,ph-M3,clk2-M3-ph,M3-M2,M2-M1,M1]; td = [td td];
%             topA = cat(3,Az,Ap1,Ap2,Ap3,Ap3,Ap2,Ap1,Az,     Az,An1,An2,An3,An3,An2,An1,Az); 
%             topB = cat(3,-B,-B,-B,-B,B,B,B,B,               B,B,B,B,-B,-B,-B,-B); 
%             Vrs = [0,1,2,3,3,2,1,0,                         0,-1,-2,-3,-3,-2,-1,0];
%             Us = [-1,-1,-1,-1,1,1,1,1,                      1,1,1,1,-1,-1,-1,-1]; 
%             inPair0 = [1,2,3,5,6,7; 2,3,4,6,7,8];   
%                 inPair02 = inPair0 + 8; inPair0 = [inPair0 inPair02]';     
%             inPolD = [1,1,1,0,0,0,                          1,1,1,0,0,0]'; 
%             inPolP = [0,0,0,0,0,0,                          0,0,0,0,0,0]';
%         end
%     else
%         error('Interval Definitions: unexpected if/else evaluation!!'); 
%     end
%     td = td*Pclk150; te = cumsum(td); te2 = [0 te]; 
% end
% % --------------------------------------------------------------







