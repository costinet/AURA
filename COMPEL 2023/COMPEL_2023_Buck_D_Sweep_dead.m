function [stick] = COMPEL_2023_Buck_D_Sweep_dead(X,Coss_adj,Ron_adj,L_adj,RL_adj)
% This function completes a steady state simulation for a  Buck
% converter
% X is a vector that contains the transistor selection, inductor selection, and frequency selection
% Coss_adj,Ron_adj,L_adj,RL_adj are the amount that the discrete components in X are perturbed for optimization
% stick is the power loss of the solved steady-state of the covnerter



% Load in the databases
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


%% Load test circuit

%modelfile = 'Buck_Vout_D.net';
modelfile = topDB.topology(6).NetList;






%% Load for Buck Converter

% Take input variables
adjust = [ones(1,2), 1, 1e-6] ;
X = X./adjust;

FETs = X(1:2);
Inductor_Chosen = X(end-1);
fs = X(end);



% FET Selection
ron = [];
Coss = [];
FET_w = [];
FET_l = [];



for i = 1:length(FETs)
    ron(i)= transDB(FETs(i)).ron.typ*1e-3;
    Coss(i) = transDB(FETs(i)).Coss.typ*1e-12;
end

% Inductor Selection
RL = indDB(Inductor_Chosen).Rdc.typ*1e-3;
L = indDB(Inductor_Chosen).L.typ*1e-6;

% Set switching interval
ON = 1;
OFF = 0;
swvec = [
    ON  OFF OFF  OFF  % Q1 ON
    OFF OFF OFF  OFF
    OFF ON  OFF  OFF  % Q2 ON
    OFF OFF OFF  OFF
    ];   % Dead

% Component Variables
Vg = 20;
Vbat1 = 4;
M = (Vbat1)/(Vg);
Vout = (Vbat1);
Pout = 25;
Dboost = 0.15;
Dbuck = M*(1-Dboost);
Vout = M*Vg;
Iout = Pout/Vout;

Rg = 1e-3;
Rbatt = 5e-3;
Rcap = 2e-3;
Cout = 4.7e-6;
RL = indDB(Inductor_Chosen).Rdc.typ*1e-3;
L = indDB(Inductor_Chosen).L.typ*1e-6;
Cfly = 23.5e-6;
CflyESR = 2.2e-3;

Vb1 = Vbat1;
Rb = Rbatt;

Ts = 1/fs;
dt = 5e-9;

Numerical_Components = {
    'C1' Cout
    'L1' L+L_adj
    'M1_C' Coss(1)+Coss_adj(1)
    'M2_C' Coss(2)+Coss_adj(2)
    'D1_C' Coss(1)+Coss_adj(1)
    'D2_C' Coss(2)+Coss_adj(2)
    'R1' RL + RL_adj
    'R2' Rbatt
    'R3' Rg
    'R4' Rcap
    'M1_R' ron(1)+Ron_adj(1)
    'M2_R' ron(2)+Ron_adj(2)
    'D1_R' ron(1)+Ron_adj(1)
    'D2_R' ron(2)+Ron_adj(2)

    };

% The inital guess of time intervals % The inital guess of time intervals
% Assigned later dynamically
ts_Baxter = [0.70  0.3];
ts = (ts_Baxter./sum(ts_Baxter)).*Ts;




% List all of the numerical components in the netlist file for all
% FETs you must use the syntax used below:

% List out all char variables in the
Switch_Resistors = {
    'M1_R'
    'M2_R'
    'D1_R'
    'D2_R'
    };

% List of the switch sequency. Organized by: the FETs (column) vs time
% interval (rows) matching Switch_Resistors and ts respectivly
Switch_Names = {
    'M1'
    'M2'
    'D1'
    'D2'
    };



ON = 1;
OFF = 0;

% List all the resistances of the diodes or FETS when they are on or
% off
SW_OFF = ones(1,4).*10000000;

%SW_ON = [M1_R,M2_R,etc...]
SW_ON = [ron(1)+Ron_adj(1) ron(2)+Ron_adj(2)  ron(1)+Ron_adj(1) ron(2)+Ron_adj(2) ];
SW = [SW_OFF;SW_ON;SW_ON];


transDB(FETs(1)).Vf.typ;

Diode_Forward_Voltage = [0 0 transDB(FETs(1)).Vf.typ transDB(FETs(2)).Vf.typ]';
u = [Vg Vbat1 0 0 transDB(FETs(1)).Vf.typ transDB(FETs(2)).Vf.typ]';
Order = [1 2 3 4];


% Initalize Variables for Modeling
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

u = u';

% Set up classes
sim = SMPSim();
conv = sim.converter;
top = sim.topology;

% Set to run numerical parsesr

% Fill in classes

conv.u = u;
top.order = Order;
top.Element_Properties = Numerical_Components;
top.Switch_Resistors = Switch_Resistors;
top.Switch_Resistor_Values = SW;
top.Switch_Sequence = swvec;
top.Fwd_Voltage = Diode_Forward_Voltage;
top.Switch_Names = Switch_Names;
%conv.ts = ts;


% Load converter into topology class
top.loadCircuit(modelfile,swvec,1);

Pout = 0;
d = M+0.01;

%% Analyze circuit
% Loop through increasing deadtime each time until Pout ~100W operating
% point
try
    while Pout < 100
        % Set duty cycle
        d = d+0.001;
        dt = 2e-9/Ts;
        ds = [d-dt  dt  1-d-dt  dt];
        ts = ds/sum(ds)*Ts;
        
        OLVin = 16/d;
        % MaxVin = OLVin+1;
        MaxVin = OLVin+2.5;
       
        Vinrange = Vg;
        % IF needed for Switched capacitor converter can change Vg to
        % affect output power for buck converter keep constant
        for Vin = Vinrange
            u(VgPOS) = Vin;
            sim.u = u;
            
            
            conv.setSwitchingPattern(1:size(swvec,1), ts)
            
            Xss = sim.steadyState;
            
            if(debug)
                sim.plotAllStates(1);
            end
            
            % Correct for diode errors
            sim.findValidSteadyState;
            Xss = sim.steadyState;
            
            [ avgXs, avgYs ] = sim.ssAvgs(Xss);
            
            % Find output and input current avg
            Ib1 = avgYs(3);
            I1 = -avgYs(1);
            
            % Find the eff, power loss and output power
            eta = (Ib1*Vb1)/(Vin*-I1);
            Ploss = -(Ib1*Vb1) + -(Vin*I1);
            Pout = (Ib1*Vb1);
            
            % Record data as needed
            etaSim = [etaSim; eta];
            PlossSim = [PlossSim; Ploss];
            PoutSim = [PoutSim; Pout];
            %             PossSim = [PossSim; Poss];
            conditions = [conditions; d, i, Vin];
            
            % If overshot the first iteration then reduce duty cycle and
            % redo
            if Pout > 100 && d == M+0.02
                d = M;
            end
            
        end
    end
    
catch ME
    J = 155465456;
end

% Set outputs
stick = Ploss;

return

end

