function [stick] = COMPEL_2023_5Phase_BuckHybridBuck_D_Sweep_GAM(X)



    load('D:\GitHub\AURA\AURAdb\databases\@transistorDB\transistors.mat')
    transDB = obj;
    clear obj

    load('D:\GitHub\AURA\AURAdb\databases\@inductorDB\inductors.mat')
    indDB = obj;
    clear obj

    load('D:\GitHub\AURA\AURAdb\@topologyDB\topology.mat')
    topDB = obj;
    clear obj


niter = 0;
debug = 0;
debug2 = 0;

if ~debug
    w = warning ('off','all');
else
    w = warning ('on','all');
end

% sdir = mfilename('fullpath');
% sdir = sdir(1:find(sdir=='\',1,'last')-1);
% addpath(sdir);

%% Load test circuit

modelfile = '5PhaseBuckHybridBuckIout.net';

%% Load for PLECS
% find_system(modelfile,'SearchDepth',1, 'IncludeCommented', 'on')
%open_system(modelfile,'loadonly');
%circuitPath = [modelfile '/' PLECsModel];
%set_param(circuitPath,'Commented','on');
%simout = sim(modelfile,eps);

%for i = 1:length(simout.properties)
%    assignin('base',simout.properties{i},eval(['simout.' simout.properties{i}]));
%end

%set_param(circuitPath,'Commented','off');




%% Load for Buck Converter

%%{
adjust = [ones(1,2), 1, 0.1, 1e-6] ;
X = X./adjust;

FETs = X(1:2);
Inductor_Chosen = X(end-2);
Iout = X(end-1);
fs = X(end);

% Number of FETs to parallel through
ParallelNumberM = 20;

% Number of Inductors to parallel through
ParallelNumberL = 20;

% FET Selection
ron = [];
Coss = [];
FET_w = [];
FET_l = [];

parallel_fix = (FETs>ParallelNumberM);
FETs(FETs>ParallelNumberM) = FETs(FETs>ParallelNumberM)-ParallelNumberM;
parallel_fix = parallel_fix + (FETs>ParallelNumberM);
FETs(FETs>ParallelNumberM) = FETs(FETs>ParallelNumberM)-ParallelNumberM;
parallel_fix = parallel_fix + (FETs>ParallelNumberM);
FETs(FETs>ParallelNumberM) = FETs(FETs>ParallelNumberM)-ParallelNumberM;

parallel_fix = parallel_fix + ones(size(parallel_fix));



for i = 1:length(FETs)
    ron(i)= (1/parallel_fix(i))*transDB(FETs(i)).ron.typ*1e-3;
    Coss(i) = (parallel_fix(i))*transDB(FETs(i)).Coss.typ*1e-12;
    Area_FET(i) = (transDB(FETs(i)).width.typ) * (transDB(FETs(i)).length.typ) * 5 * parallel_fix(i);
end


parallel_fixL = (Inductor_Chosen>ParallelNumberL);
Inductor_Chosen(Inductor_Chosen>ParallelNumberL) = Inductor_Chosen(Inductor_Chosen>ParallelNumberL)-ParallelNumberL;
parallel_fixL = parallel_fixL + (Inductor_Chosen>ParallelNumberL);
Inductor_Chosen(Inductor_Chosen>ParallelNumberL) = Inductor_Chosen(Inductor_Chosen>ParallelNumberL)-ParallelNumberL;
parallel_fixL = parallel_fixL + (Inductor_Chosen>ParallelNumberL);
Inductor_Chosen(Inductor_Chosen>ParallelNumberL) = Inductor_Chosen(Inductor_Chosen>ParallelNumberL)-ParallelNumberL;

parallel_fixL = parallel_fixL + ones(size(parallel_fixL));

RL = (1/parallel_fixL)*indDB(Inductor_Chosen).Rdc.typ*1e-3;
L = (1/parallel_fixL)*indDB(Inductor_Chosen).L.typ*1e-6;
Area_L = (indDB(Inductor_Chosen).width.typ) * (indDB(Inductor_Chosen).length.typ) * 5 * parallel_fixL;

% Set switching interval
ON = 1;
OFF = 0;
swvec = [
    1 0 1 0 1 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 ;
    0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 ; 
    0 1 0 1 0 1 0 1 0 1 0 0 0 0 0 0 0 0 0 0 ;
    0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 ; 
    ];   % Dead


Vg = 48;
Vbat1 = 1;
M = (Vbat1)/(Vg);
Vout = (Vbat1);
Pout = 25;
Dboost = 0.15;
Dbuck = M*(1-Dboost);
Vout = M*Vg;
% Iout = 10; Set by function

Rg = 1e-3;
Rbatt = 5e-3;
Rcap = 4e-3;
Cout = 9.4e-6;
Cfly = 9.4e-6;
CflyESR = 4e-3;

Vb1 = Vbat1;
Rb = Rbatt;

Ts = 1/fs;
dt = 5e-9;

Coss_adj = [0 0];
Ron_adj = [0 0];
L_adj = [0];
RL_adj = [0];

Numerical_Components = {
    'C6' Cout
    'L1' L+L_adj
    'L2' L+L_adj
    'L3' L+L_adj
    'L4' L+L_adj
    'L5' L+L_adj
    'R4' RL + RL_adj
    'R6' RL + RL_adj
    'R8' RL + RL_adj
    'R7' RL + RL_adj
    'R9' RL + RL_adj
    'C1' Cfly
    'C2' Cfly
    'C3' Cfly
    'C4' Cfly
    'R2' CflyESR
    'R3' CflyESR
    'R10' CflyESR
    'R11' CflyESR
    'R5' Rbatt
    'R1' Rg
    'M1_C' Coss(1) + Coss_adj(1)
    'M2_C' Coss(1) + Coss_adj(1)
    'M3_C' Coss(1) + Coss_adj(1)
    'M4_C' Coss(1) + Coss_adj(1)
    'M5_C' Coss(1) + Coss_adj(1)
    'M6_C' Coss(2) + Coss_adj(2)
    'M7_C' Coss(2) + Coss_adj(2)
    'M8_C' Coss(2) + Coss_adj(2)
    'M9_C' Coss(2) + Coss_adj(2)
    'M10_C' Coss(2) + Coss_adj(2)
    'D1_C' Coss(1) + Coss_adj(1)
    'D2_C' Coss(1) + Coss_adj(1)
    'D3_C' Coss(1) + Coss_adj(1)
    'D4_C' Coss(1) + Coss_adj(1)
    'D5_C' Coss(1) + Coss_adj(1)
    'D6_C' Coss(2) + Coss_adj(2)
    'D7_C' Coss(2) + Coss_adj(2)
    'D8_C' Coss(2) + Coss_adj(2)
    'D9_C' Coss(2) + Coss_adj(2)
    'D10_C' Coss(2) + Coss_adj(2)
    'M1_R' ron(1) + Ron_adj(1)
    'M2_R' ron(1) + Ron_adj(1)
    'M3_R' ron(1) + Ron_adj(1)
    'M4_R' ron(1) + Ron_adj(1)
    'M5_R' ron(1) + Ron_adj(1)
    'M6_R' ron(2) + Ron_adj(2)
    'M7_R' ron(2) + Ron_adj(2)
    'M8_R' ron(2) + Ron_adj(2)
    'M9_R' ron(2) + Ron_adj(2)
    'M10_R' ron(2) + Ron_adj(2)
    'D1_R' ron(1) + Ron_adj(1)
    'D2_R' ron(1) + Ron_adj(1)
    'D3_R' ron(1) + Ron_adj(1)
    'D4_R' ron(1) + Ron_adj(1)
    'D5_R' ron(1) + Ron_adj(1)
    'D6_R' ron(2) + Ron_adj(2)
    'D7_R' ron(2) + Ron_adj(2)
    'D8_R' ron(2) + Ron_adj(2)
    'D9_R' ron(2) + Ron_adj(2)
    'D10_R' ron(2) + Ron_adj(2)
    };

ts_Baxter = [0.2  0.8  0.2  0.8 ];
ts = (ts_Baxter./sum(ts_Baxter)).*Ts;

% The inital guess of time intervals % The inital guess of time intervals
% Assigned later dynamically


% List all of the numerical components in the netlist file for all
% FETs you must use the syntax used below:


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
    'D1_R'
    'D2_R'
    'D3_R'
    'D4_R'
    'D5_R'
    'D6_R'
    'D7_R'
    'D8_R'
    'D9_R'
    'D10_R'
    };
% List of the switch sequency. Organized by: the FETs (column) vs time
% interval (rows) matching Switch_Resistors and ts respectivly

Switch_Names = {
    'M1'
    'M2'
    'M3'
    'M4'
    'M5'
    'M6'
    'M7'
    'M8'
    'M9'
    'M10'
    'D1'
    'D2'
    'D3'
    'D4'
    'D5'
    'D6'
    'D7'
    'D8'
    'D9'
    'D10'
    };




ON = 1;
OFF = 0;

% List all the resistances of the diodes or FETS when they are on or
% off
SW_OFF = ones(1,20).*10000000;

%SW_ON = [M1_R,M2_R,etc...]
SW_ON = [ron(1) + Ron_adj(1)
    ron(1) + Ron_adj(1)
    ron(1) + Ron_adj(1)
    ron(1) + Ron_adj(1)
    ron(1) + Ron_adj(1)
    ron(2) + Ron_adj(2)
    ron(2) + Ron_adj(2)
    ron(2) + Ron_adj(2)
    ron(2) + Ron_adj(2)
    ron(2) + Ron_adj(2)
    ron(1) + Ron_adj(1)
    ron(1) + Ron_adj(1)
    ron(1) + Ron_adj(1)
    ron(1) + Ron_adj(1)
    ron(1) + Ron_adj(1)
    ron(2) + Ron_adj(2)
    ron(2) + Ron_adj(2)
    ron(2) + Ron_adj(2)
    ron(2) + Ron_adj(2)
    ron(2) + Ron_adj(2)
    ]';

SW = [SW_OFF;SW_ON;SW_ON];




Diode_Forward_Voltage = [0
    0
    0
    0
    0
    0
    0
    0
    0
    0
    transDB(FETs(1)).Vf.typ
    transDB(FETs(1)).Vf.typ
    transDB(FETs(1)).Vf.typ
    transDB(FETs(1)).Vf.typ
    transDB(FETs(1)).Vf.typ
    transDB(FETs(2)).Vf.typ
    transDB(FETs(2)).Vf.typ
    transDB(FETs(2)).Vf.typ
    transDB(FETs(2)).Vf.typ
    transDB(FETs(2)).Vf.typ
    ];

u = [Vg
    0
    0
    0
    0
    0
    0
    0
    0
    0
    0
    transDB(FETs(1)).Vf.typ
    transDB(FETs(1)).Vf.typ
    transDB(FETs(1)).Vf.typ
    transDB(FETs(1)).Vf.typ
    transDB(FETs(1)).Vf.typ
    transDB(FETs(2)).Vf.typ
    transDB(FETs(2)).Vf.typ
    transDB(FETs(2)).Vf.typ
    transDB(FETs(2)).Vf.typ
    transDB(FETs(2)).Vf.typ
    Iout
    ]';

Order = [1 2 3 4];
%}
etaSim = [];
PlossSim = [];
PoutSim = [];
conditions = [];
PossSim = [];
drange = linspace(0.8, .99, 12);
Vscloc1 = 15;
Vscloc2 = 18;
Illoc = 1;
Ib1loc = 2;
Ib2loc = 3;
VgPOS = 1;
i = 1;



sim = SMPSim();
conv = sim.converter;
top = sim.topology;

% Set to run numerical parsesr

%conv.ts = ts;
conv.u = u;
top.order = Order;
top.Element_Properties = Numerical_Components;
top.Switch_Resistors = Switch_Resistors;
top.Switch_Resistor_Values = SW;
top.Switch_Sequence = swvec;
top.Fwd_Voltage = Diode_Forward_Voltage;
top.Switch_Names = Switch_Names;
%



top.loadCircuit(modelfile,swvec,1);

Pout = 0;
d = 0.21;
starting_d = d;
%% Analyze circuit
try
    while Pout < Iout
       
        ds = [d 1-d d 1-d];
        ts = ds/sum(ds)*Ts;
        
        OLVin = 16/d;
        % MaxVin = OLVin+1;
        MaxVin = OLVin+2.5;
       
        Vinrange = Vg;
        for Vin = Vinrange
            u(VgPOS) = Vin;
            sim.u = u;
            
            
            conv.setSwitchingPattern(1:size(swvec,1), ts)
            
            Xss = sim.steadyState;
         if (debug)
                sim.plotAllStates(1);
        end
            
            % ssOrder = plecs('get', circuitPath, 'StateSpaceOrder');
            % outputs = ssOrder.Outputs;
            
            outputs = top.outputLabels;
            
            %% finalRun
            % once everything seems to be error-free based on discrete time points,
            % goes through once more with eigenvalue-based spacing to make sure no
            % inter-sample violations are occuring.
            finalRun = 0;
            
            %% Symmetry check
            % May be useful but not doing anything with it yet.  Can identify that DAB,
            % etc. exhibit half-cycle symmetry
            % TF = conv.checkForSymmetry;
            
            %%
            %[Xf,ts,swinds] = timeSteppingPeriod(sim);

            
           % DMC_SteadyState(sim,conv,top);
            Xss = sim.steadyState;
            
            [ avgXs, avgYs ] = sim.ssAvgs(Xss);
            
            Ib1 = Iout;
            Vb1 = avgYs(48)+avgYs(18)*Rbatt;
           % Ib2 = (avgYs(25)+avgYs(8)*Rcap-Vbat2)/Rbatt;
            I1 = -avgYs(1);
            
            
            eta = (Ib1*Vb1)/(Vin*-I1);
            Ploss = -(Ib1*Vb1) + -(Vin*I1);
            Pout = (Ib1*Vb1);
            
          
            etaSim = [etaSim; eta];
            PlossSim = [PlossSim; Ploss];
            PoutSim = [PoutSim; Pout];
            %             PossSim = [PossSim; Poss];
            conditions = [conditions; d, i, Vin];
            
            if Pout > Iout && d == starting_d
                d = 0.20;
                Pout = 0;
            end
            
            
            if Pout <Iout-0.2
                d = d+0.00005;
            else
                d = d+0.00001;
            end
            
            
        end
    end
    
catch ME
    J = 155465456;
end



MMperA = (sum(Area_FET)+Area_L)/Iout;
ONEminusEff = (1-eta);
stick = [MMperA ONEminusEff];

% Saving each iteration
load('COMPEL_2023_5PhaseBuckHybridBuck_D_Sweep_GA_saved.mat','sav_x','sav_eff','sav_Ploss','sav_Pout','sav_conditions','sav_area','sav_stick');

sav_x = [sav_x; X];
sav_eff = [sav_eff; eta];
sav_Ploss = [sav_Ploss; Ploss];
sav_Pout = [sav_Pout; Pout];
sav_conditions = [sav_conditions; d,Vb1,Vin];
sav_area = [sav_area ; Area_FET Area_L];
sav_stick = [sav_stick; stick];


save('COMPEL_2023_5PhaseBuckHybridBuck_D_Sweep_GA_saved.mat','sav_x','sav_eff','sav_Ploss','sav_Pout','sav_conditions','sav_area','sav_stick');

return



end

