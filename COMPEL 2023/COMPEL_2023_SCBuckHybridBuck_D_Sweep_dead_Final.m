function [stick] = COMPEL_2023_SCBuckHybridBuck_D_Sweep_dead_Final(X,Coss_adj,Ron_adj,L_adj,RL_adj)

try
    load('D:\GitHub\AURA\AURAdb\databases\@transistorDB\transistors.mat')
    transDB = obj;
    clear obj
    
    load('D:\GitHub\AURA\AURAdb\databases\@inductorDB\inductorsFull.mat')
    indDB = obj;
    clear obj
    
    load('D:\GitHub\AURA\AURAdb\@topologyDB\topology.mat')
    topDB = obj;
    clear obj

    load('D:\GitHub\AURA\AURAdb\databases\@capacitorDB\capacitors.mat')
    capDB = obj;
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

modelfile = 'SCBuckHybridBuckIout.net';


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




%% Load Converter

%%{
adjust = [ones(1,5), 1, 1e-6  ones(1,5) ones(1,1) ones(1,2) 0.1] ;
X = X./adjust;

FETs = X(1:5);
Inductor_Chosen = X(6);
fs = X(7);
PFETs = X(8:12);
PINDs = X(13);
PCfly = X(14);
PCbuck = X(15);
Iout = X(16);



% Number of FETs to parallel through
ParallelNumberM = 20;

% Number of Inductors to parallel through
ParallelNumberL = 16;

% FET Selection
ron = [];
Coss = [];
FET_w = [];
FET_l = [];

for i = 1:length(FETs)
    %{
    if i == 1 && FETs(i)==6
        stick = [100  1];
        return
    end
    if i == 3  && FETs(i)>5
        stick = [100  1];
        return
    end
    if i == 4  && FETs(i)>5
        stick = [100  1];
        return
    end
    if i == 5  && FETs(i)>5
        stick = [100  1];
        return
    end
    %}
    ron(i)= (1/PFETs(i))*transDB(FETs(i)).ron.typ*1e-3;
    Coss(i) = (PFETs(i))*transDB(FETs(i)).Coss.typ*1e-12;
    Area_FET(i) = (transDB(FETs(i)).width.typ) * (transDB(FETs(i)).length.typ) * 3 * PFETs(i);
end

RL = (1/PINDs)*indDB(Inductor_Chosen).Rdc.typ*1e-3;
L = (1/PINDs)*indDB(Inductor_Chosen).L.typ*1e-6;
Area_L = (indDB(Inductor_Chosen).width.typ) * (indDB(Inductor_Chosen).length.typ) * 3 * PINDs;




% Set switching interval
ON = 1;
OFF = 0;
swvec = [
    [1 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; %dead
    0 0 0 1 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; %dead
    0 1 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; %dead
    0 0 0 1 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; %dead
    ] + ...
    [0 0 0 0 1 0 0 0 0 1 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 1 0 0 0 0 1 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; %dead
    0 0 0 0 1 0 0 0 0 1 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; %dead
    0 0 0 0 0 1 1 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 1 1 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; %dead
    0 0 0 0 0 1 1 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; %dead
    ]
    ];


Vg = 12;
Vbat1 = 1;
M = (Vbat1)/(Vg);
Vout = (Vbat1);
Pout = 25;
Dboost = 0.15;
Dbuck = M*(1-Dboost);
Vout = M*Vg;
%Iout = 10; % Now set by function

%{
Rg = 1e-3;
Rbatt = 1e-4;
Rcap = 0.1e-3;
Cout = 10*9.4e-6;
Cfly = 10*9.4e-6;
CflyESR = 0.1e-3;
%}
%{
Rg = 1e-3;
Rbatt = 5e-3;
Rcap = 4e-3;
Cout = 10*9.4e-6;
Cfly = 10*9.4e-6;
CflyESR = 4e-3;
%}

Rg = 1e-3;
Rbatt = 1e-4;
Rcap = 0.1e-3;
Cout = 200e-6;

Cfly = PCfly*capDB(16).C.typ*1e-6;
CflyESR = (1/PCfly)*capDB(16).ESR.typ*1e-3;
Area_C = PCfly*capDB(16).width.typ*capDB(16).length.typ*3;

Cbuck = PCbuck*capDB(9).C.typ*1e-6;
CbuckESR = (1/PCbuck)*capDB(9).ESR.typ*1e-3;
Area_C = Area_C + PCbuck*capDB(9).width.typ*capDB(9).length.typ*2;


Vb1 = Vbat1;
Rb = Rbatt;

Ts = 1/fs;
dt = 5e-9;
%Coss_adj = [0 0 0 0 0];
%Ron_adj = [0 0 0 0 0];
%L_adj = [0];
%RL_adj = [0];

Numerical_Components = {
    'C6' Cout
    'L1' L+L_adj
    'L2' L+L_adj
    'L3' L+L_adj
    'R4' RL + RL_adj
    'R6' RL + RL_adj
    'R8' RL + RL_adj
    'C1' Cbuck
    'C2' Cbuck
    'R2' CbuckESR
    'R3' CbuckESR
    'C3' Cfly
    'C4' Cfly
    'C5' Cfly
    'R5' CflyESR
    'R9' CflyESR
    'R7' CflyESR
    'R1' Rg
    'R10' Rcap
    'M1_C' Coss(1) + Coss_adj(1)
    'M2_C' Coss(1) + Coss_adj(1)
    'M3_C' Coss(1) + Coss_adj(1)
    'M4_C' Coss(2) + Coss_adj(2)
    'M5_C' Coss(3) + Coss_adj(3)
    'M6_C' Coss(4) + Coss_adj(4)
    'M7_C' Coss(5) + Coss_adj(5)
    'M8_C' Coss(2) + Coss_adj(2)
    'M9_C' Coss(2) + Coss_adj(2)
    'M10_C' Coss(3) + Coss_adj(3)
    'M11_C' Coss(4) + Coss_adj(4)
    'M12_C' Coss(5) + Coss_adj(5)
    'M13_C' Coss(3) + Coss_adj(3)
    'M14_C' Coss(4) + Coss_adj(4)
    'M15_C' Coss(5) + Coss_adj(5)
    'D1_C' Coss(1) + Coss_adj(1)
    'D2_C' Coss(1) + Coss_adj(1)
    'D3_C' Coss(1) + Coss_adj(1)
    'D4_C' Coss(2) + Coss_adj(2)
    'D5_C' Coss(3) + Coss_adj(3)
    'D6_C' Coss(4) + Coss_adj(4)
    'D7_C' Coss(5) + Coss_adj(5)
    'D8_C' Coss(2) + Coss_adj(2)
    'D9_C' Coss(2) + Coss_adj(2)
    'D10_C' Coss(3) + Coss_adj(3)
    'D11_C' Coss(4) + Coss_adj(4)
    'D12_C' Coss(5) + Coss_adj(5)
    'D13_C' Coss(3) + Coss_adj(3)
    'D14_C' Coss(4) + Coss_adj(4)
    'D15_C' Coss(5) + Coss_adj(5)
    'M1_R' ron(1) + Ron_adj(1)
    'M2_R' ron(1) + Ron_adj(1)
    'M3_R' ron(1) + Ron_adj(1)
    'M4_R' ron(2) + Ron_adj(2)
    'M5_R' ron(3) + Ron_adj(3)
    'M6_R' ron(4) + Ron_adj(4)
    'M7_R' ron(5) + Ron_adj(5)
    'M8_R' ron(2) + Ron_adj(2)
    'M9_R' ron(2) + Ron_adj(2)
    'M10_R' ron(3) + Ron_adj(3)
    'M11_R' ron(4) + Ron_adj(4)
    'M12_R' ron(5) + Ron_adj(5)
    'M13_R' ron(3) + Ron_adj(3)
    'M14_R' ron(4) + Ron_adj(4)
    'M15_R' ron(5) + Ron_adj(5)
    'D1_R' ron(1) + Ron_adj(1)
    'D2_R' ron(1) + Ron_adj(1)
    'D3_R' ron(1) + Ron_adj(1)
    'D4_R' ron(2) + Ron_adj(2)
    'D5_R' ron(3) + Ron_adj(3)
    'D6_R' ron(4) + Ron_adj(4)
    'D7_R' ron(5) + Ron_adj(5)
    'D8_R' ron(2) + Ron_adj(2)
    'D9_R' ron(2) + Ron_adj(2)
    'D10_R' ron(3) + Ron_adj(3)
    'D11_R' ron(4) + Ron_adj(4)
    'D12_R' ron(5) + Ron_adj(5)
    'D13_R' ron(3) + Ron_adj(3)
    'D14_R' ron(4) + Ron_adj(4)
    'D15_R' ron(5) + Ron_adj(5)
    };

ts_Baxter = [0.2  0.8  0.2  0.8 0.2  0.8  0.2  0.8 ];
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
    'M11_R'
    'M12_R'
    'M13_R'
    'M14_R'
    'M15_R'
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
    'D11_R'
    'D12_R'
    'D13_R'
    'D14_R'
    'D15_R'
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
    'M11'
    'M12'
    'M13'
    'M14'
    'M15'
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
    'D11'
    'D12'
    'D13'
    'D14'
    'D15'
    };




ON = 1;
OFF = 0;

% List all the resistances of the diodes or FETS when they are on or
% off
SW_OFF = ones(1,30).*10000000;

%SW_ON = [M1_R,M2_R,etc...]
SW_ON = [   ron(1) + Ron_adj(1)
    ron(1) + Ron_adj(1)
    ron(1) + Ron_adj(1)
    ron(2) + Ron_adj(2)
    ron(3) + Ron_adj(3)
    ron(4) + Ron_adj(4)
    ron(5) + Ron_adj(5)
    ron(2) + Ron_adj(2)
    ron(2) + Ron_adj(2)
    ron(3) + Ron_adj(3)
    ron(4) + Ron_adj(4)
    ron(5) + Ron_adj(5)
    ron(3) + Ron_adj(3)
    ron(4) + Ron_adj(4)
    ron(5) + Ron_adj(5)
    ron(1) + Ron_adj(1)
    ron(1) + Ron_adj(1)
    ron(1) + Ron_adj(1)
    ron(2) + Ron_adj(2)
    ron(3) + Ron_adj(3)
    ron(4) + Ron_adj(4)
    ron(5) + Ron_adj(5)
    ron(2) + Ron_adj(2)
    ron(2) + Ron_adj(2)
    ron(3) + Ron_adj(3)
    ron(4) + Ron_adj(4)
    ron(5) + Ron_adj(5)
    ron(3) + Ron_adj(3)
    ron(4) + Ron_adj(4)
    ron(5) + Ron_adj(5)
]';

SW = [SW_OFF;SW_ON;SW_ON];


transDB(FETs(1)).Vf.typ;

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
    0
    0
    0
    0
    0
    transDB(FETs(1)).Vf.typ
    transDB(FETs(1)).Vf.typ
    transDB(FETs(1)).Vf.typ
    transDB(FETs(2)).Vf.typ
    transDB(FETs(3)).Vf.typ
    transDB(FETs(4)).Vf.typ
    transDB(FETs(5)).Vf.typ
    transDB(FETs(2)).Vf.typ
    transDB(FETs(2)).Vf.typ
    transDB(FETs(3)).Vf.typ
    transDB(FETs(4)).Vf.typ
    transDB(FETs(5)).Vf.typ
    transDB(FETs(3)).Vf.typ
    transDB(FETs(4)).Vf.typ
    transDB(FETs(5)).Vf.typ
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
    0
    0
    0
    0
    0
    transDB(FETs(1)).Vf.typ
    transDB(FETs(1)).Vf.typ
    transDB(FETs(1)).Vf.typ
    transDB(FETs(2)).Vf.typ
    transDB(FETs(3)).Vf.typ
    transDB(FETs(4)).Vf.typ
    transDB(FETs(5)).Vf.typ
    transDB(FETs(2)).Vf.typ
    transDB(FETs(2)).Vf.typ
    transDB(FETs(3)).Vf.typ
    transDB(FETs(4)).Vf.typ
    transDB(FETs(5)).Vf.typ
    transDB(FETs(3)).Vf.typ
    transDB(FETs(4)).Vf.typ
    transDB(FETs(5)).Vf.typ
    Iout]';

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
d = 0.75;
starting_d = d;

%% Analyze circuit

    while Pout < Iout

        
        dt = 0.00000080;
        ds = [d-dt  dt  1-d-dt  dt   d-dt   dt    1-d-dt  dt];
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
            if(debug)
                sim.plotAllStates(1);
                sim.plotAllOutputs(2);
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

            
           % sim.findValidSteadyState;
            %COMPEL_2023_SteadyState(sim,conv,top)
            Xss = sim.steadyState;
            
            [ avgXs, avgYs ] = sim.ssAvgs(Xss);
            
            
            Vb1 = avgYs(78)+avgYs(39)*Rbatt;
            Ib1 = Iout;
            I1 = -avgYs(1);
            
            %% Magnetic Losses
            %%{
            alpha = 2;
            beta = 2.75;
            Kfe = 4e-6;
            number_of_turns = 1;
            le = 0.7854; % in cm 2.5mm x pi
            Ae = 0.2275; % in cm3 7mm x 5 mm x 6.5 mm
            Ae_in_m = Ae/100/100/100;
            DeltaH = number_of_turns*(max(Xss(20,:))-min(Xss(20,:)))/le;
                
            PL=Kfe*(fs)^alpha*(DeltaH)^beta; % in kW/m3

            Core_Loss = PL/1000*Ae_in_m*3;
            %}
           % Core_Loss = 0;

            eta = (Ib1*Vb1)/((Vin*-I1));
            Ploss = -(Ib1*Vb1) + -(Vin*I1);
            Pout = (Ib1*Vb1);
            
           
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

            
            if Pout > Iout && d == starting_d
                d = 0.18;
            end
            
            if Pout<Iout-0.1
                d = d+0.0005;
                if Pout<Iout-2
                    d = d+0.005;
                end
            else
                d = d+0.00005;
            end
            
        end
    end
    
catch ME
    J = 155465456;
end


MMperA = (sum(Area_FET)+Area_L+sum(Area_C))/Iout;
ONEminusEff = (1-eta);
stick = [MMperA ONEminusEff Core_Loss Pout Ploss];

% 
% % Saving each iteration
% load('COMPEL_2023_SCBuckHybridBuck_D_Sweep_GA_saved_C2.mat','sav_x','sav_eff','sav_Ploss','sav_Pout','sav_conditions','sav_area','sav_stick');
% 
% sav_x = [sav_x; X];
% sav_eff = [sav_eff; eta];
% sav_Ploss = [sav_Ploss; Ploss];
% sav_Pout = [sav_Pout; Pout];
% sav_conditions = [sav_conditions; d,Vb1,Vin];
% sav_area = [sav_area ; Area_FET Area_L Area_C];
% sav_stick = [sav_stick; stick];
% 
% 
% save('COMPEL_2023_SCBuckHybridBuck_D_Sweep_GA_saved_C2.mat','sav_x','sav_eff','sav_Ploss','sav_Pout','sav_conditions','sav_area','sav_stick');

return

end

