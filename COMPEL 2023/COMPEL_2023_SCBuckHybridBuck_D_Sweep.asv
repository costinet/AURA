function [stick,graph_values] = COMPEL_2023_SCBuckHybridBuck_D_Sweep(X,Coss_adj,Ron_adj,L_adj,RL_adj)


    load('D:\GitHub\AURA\AURAdb\databases\@transistorDB\transistors.mat')
    transDBig = obj;
    clear obj
    
    load('D:\GitHub\AURA\AURAdb\databases\@inductorDB\inductors.mat')
    indDBig = obj;
    clear obj
    
    load('D:\GitHub\AURA\AURAdb\@topologyDB\topology.mat')
    topDB = obj;
    clear obj
indDBindex = [    
   1
   0
   0
   0
   1
   0
   0
   0
   0
   0
   0
   0
   1
   1
   0
   0
   1
   0
   0
   0];


transDBindex = [    
   0
   0
   0
   0
   0
   0
   0
   0
   0
   1
   1
   1
   0
   0
   0
   0
   0
   1
   1
   1];
    

transDB = transistorDB;
transDB.add(transDBig(logical(transDBindex)));

indDB = inductorDB;
indDB.add(indDBig(logical(indDBindex)));


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

modelfile = topDB.topology(1).NetList;


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
adjust = [ones(1,5), 1, 1e-6] ;
X = X./adjust;

FETs = X(1:5);
Inductor_Chosen = X(end-1);
fs = X(end);


% FET Selection
ron = [];
Coss = [];
FET_w = [];
FET_l = [];

ron(1)= transDB(FETs(1)).ron.typ*1e-3;
Coss(1) = transDB(FETs(1)).Coss.typ*1e-12;
ron(2) = transDB(FETs(2)).ron.typ*1e-3;
Coss(2) = transDB(FETs(2)).Coss.typ*1e-12;
ron(3) = transDB(FETs(3)).ron.typ*1e-3;
Coss(3) = transDB(FETs(3)).Coss.typ*1e-12;
ron(4) = transDB(FETs(4)).ron.typ*1e-3;
Coss(4) = transDB(FETs(4)).Coss.typ*1e-12;
ron(5) = transDB(FETs(5)).ron.typ*1e-3;
Coss(5) = transDB(FETs(5)).Coss.typ*1e-12;

% Set switching interval
ON = 1;
OFF = 0;
swvec = [
    [1 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 1 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 1 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 1 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; ] + ...
    [0 0 0 0 1 0 0 0 0 1 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 1 0 0 0 0 1 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 1 1 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 1 1 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;  ]
    ];   % Dead


Vg = 48;
Vbat1 = 1;
M = (Vbat1)/(Vg);
Vout = (Vbat1);
Pout = 25;
Dboost = 0.15;
Dbuck = M*(1-Dboost);
Vout = M*Vg;
Iout = Pout/Vout;

Rg = 1e-3;
Rbatt = 5e-3;
Rcap = 4e-3;
Cout = 9.4e-6;
RL = indDB(Inductor_Chosen).Rdc.typ*1e-3;
L = indDB(Inductor_Chosen).L.typ*1e-6;
Cfly = 9.4e-6;
CflyESR = 4e-3;

Vb1 = Vbat1;
Rb = Rbatt;

Ts = 1/fs;
dt = 5e-9;

Numerical_Components = {
    'C6' Cout
    'L1' L+L_adj
    'L2' L+L_adj
    'L3' L+L_adj
    'R4' RL + RL_adj
    'R6' RL + RL_adj
    'R8' RL + RL_adj
    'C1' Cfly
    'C2' Cfly
    'R2' CflyESR
    'R3' CflyESR
    'C3' Cfly
    'C4' Cfly
    'C5' Cfly
    'R5' CflyESR
    'R9' CflyESR
    'R7' CflyESR
    'R1' Rg
    'R10' Rcap
    'R11' Rbatt
    'M1_C' Coss(1) + Coss_adj(1)
    'M2_C' Coss(2) + Coss_adj(2)
    'M3_C' Coss(2) + Coss_adj(2)
    'M4_C' Coss(1) + Coss_adj(1)
    'M5_C' Coss(3) + Coss_adj(3)
    'M6_C' Coss(4) + Coss_adj(4)
    'M7_C' Coss(5) + Coss_adj(5)
    'M8_C' Coss(1) + Coss_adj(1)
    'M9_C' Coss(1) + Coss_adj(1)
    'M10_C' Coss(3) + Coss_adj(3)
    'M11_C' Coss(4) + Coss_adj(4)
    'M12_C' Coss(5) + Coss_adj(5)
    'M13_C' Coss(3) + Coss_adj(3)
    'M14_C' Coss(4) + Coss_adj(4)
    'M15_C' Coss(5) + Coss_adj(5)
    'D1_C' Coss(1) + Coss_adj(1)
    'D2_C' Coss(2) + Coss_adj(2)
    'D3_C' Coss(2) + Coss_adj(2)
    'D4_C' Coss(1) + Coss_adj(1)
    'D5_C' Coss(3) + Coss_adj(3)
    'D6_C' Coss(4) + Coss_adj(4)
    'D7_C' Coss(5) + Coss_adj(5)
    'D8_C' Coss(1) + Coss_adj(1)
    'D9_C' Coss(1) + Coss_adj(1)
    'D10_C' Coss(3) + Coss_adj(3)
    'D11_C' Coss(4) + Coss_adj(4)
    'D12_C' Coss(5) + Coss_adj(5)
    'D13_C' Coss(3) + Coss_adj(3)
    'D14_C' Coss(4) + Coss_adj(4)
    'D15_C' Coss(5) + Coss_adj(5)
    'M1_R' ron(1) + Ron_adj(1)
    'M2_R' ron(2) + Ron_adj(2)
    'M3_R' ron(2) + Ron_adj(2)
    'M4_R' ron(1) + Ron_adj(1)
    'M5_R' ron(3) + Ron_adj(3)
    'M6_R' ron(4) + Ron_adj(4)
    'M7_R' ron(5) + Ron_adj(5)
    'M8_R' ron(1) + Ron_adj(1)
    'M9_R' ron(1) + Ron_adj(1)
    'M10_R' ron(3) + Ron_adj(3)
    'M11_R' ron(4) + Ron_adj(4)
    'M12_R' ron(5) + Ron_adj(5)
    'M13_R' ron(3) + Ron_adj(3)
    'M14_R' ron(4) + Ron_adj(4)
    'M15_R' ron(5) + Ron_adj(5)
    'D1_R' ron(1) + Ron_adj(1)
    'D2_R' ron(2) + Ron_adj(2)
    'D3_R' ron(2) + Ron_adj(2)
    'D4_R' ron(1) + Ron_adj(1)
    'D5_R' ron(3) + Ron_adj(3)
    'D6_R' ron(4) + Ron_adj(4)
    'D7_R' ron(5) + Ron_adj(5)
    'D8_R' ron(1) + Ron_adj(1)
    'D9_R' ron(1) + Ron_adj(1)
    'D10_R' ron(3) + Ron_adj(3)
    'D11_R' ron(4) + Ron_adj(4)
    'D12_R' ron(5) + Ron_adj(5)
    'D13_R' ron(3) + Ron_adj(3)
    'D14_R' ron(4) + Ron_adj(4)
    'D15_R' ron(5) + Ron_adj(5)
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
    ron(2) + Ron_adj(2)
    ron(2) + Ron_adj(2)
    ron(1) + Ron_adj(1)
    ron(3) + Ron_adj(3)
    ron(4) + Ron_adj(4)
    ron(5) + Ron_adj(5)
    ron(1) + Ron_adj(1)
    ron(1) + Ron_adj(1)
    ron(3) + Ron_adj(3)
    ron(4) + Ron_adj(4)
    ron(5) + Ron_adj(5)
    ron(3) + Ron_adj(3)
    ron(4) + Ron_adj(4)
    ron(5) + Ron_adj(5)
    ron(1) + Ron_adj(1)
    ron(2) + Ron_adj(2)
    ron(2) + Ron_adj(2)
    ron(1) + Ron_adj(1)
    ron(3) + Ron_adj(3)
    ron(4) + Ron_adj(4)
    ron(5) + Ron_adj(5)
    ron(1) + Ron_adj(1)
    ron(1) + Ron_adj(1)
    ron(3) + Ron_adj(3)
    ron(4) + Ron_adj(4)
    ron(5) + Ron_adj(5)
    ron(3) + Ron_adj(3)
    ron(4) + Ron_adj(4)
    ron(5) + Ron_adj(5)]';

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
    transDB(FETs(2)).Vf.typ
    transDB(FETs(2)).Vf.typ
    transDB(FETs(1)).Vf.typ
    transDB(FETs(3)).Vf.typ
    transDB(FETs(4)).Vf.typ
    transDB(FETs(5)).Vf.typ
    transDB(FETs(1)).Vf.typ
    transDB(FETs(1)).Vf.typ
    transDB(FETs(3)).Vf.typ
    transDB(FETs(4)).Vf.typ
    transDB(FETs(5)).Vf.typ
    transDB(FETs(3)).Vf.typ
    transDB(FETs(4)).Vf.typ
    transDB(FETs(5)).Vf.typ
    ];

u = [Vg 
    Vbat1
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
    transDB(FETs(2)).Vf.typ
    transDB(FETs(2)).Vf.typ
    transDB(FETs(1)).Vf.typ
    transDB(FETs(3)).Vf.typ
    transDB(FETs(4)).Vf.typ
    transDB(FETs(5)).Vf.typ
    transDB(FETs(1)).Vf.typ
    transDB(FETs(1)).Vf.typ
    transDB(FETs(3)).Vf.typ
    transDB(FETs(4)).Vf.typ
    transDB(FETs(5)).Vf.typ
    transDB(FETs(3)).Vf.typ
    transDB(FETs(4)).Vf.typ
    transDB(FETs(5)).Vf.typ]';

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
d = 0.31-0.005;

%% Analyze circuit
try
    while Pout < 100
        if Pout<90
            d = d+0.005;
        else
            d = d+0.0005;
        end
        
        % dt = 0.0025;
        ds = [d  1-d  d 1-d];
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

            
            DMC_SteadyState(sim,conv,top);
            Xss = sim.steadyState;
            
            [ avgXs, avgYs ] = sim.ssAvgs(Xss);
            
            
            Ib1 = (avgYs(78)+avgYs(39)*Rcap-Vbat1)/Rbatt;
           
            I1 = -avgYs(1);
            
            
            eta = (Ib1*Vb1)/(Vin*-I1);
            Ploss = -(Ib1*Vb1) + -(Vin*I1);
            Pout = (Ib1*Vb1);
            
            
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
            
            if Pout > 100 && d == 0.31
                d = 0.29;
            end
            
        end
    end
    
catch ME
    J = 155465456;
end

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mineff = .3;
etaSim(etaSim<mineff) = mineff;
etaSim(etaSim>1) = mineff;

locs = etaSim > mineff;

VgSim = conditions(locs,end);
F = scatteredInterpolant(VgSim, PoutSim(locs), PlossSim(locs), 'linear','none');

figure(106)

Vgrange = 14:.25:22;
PoutRange = 0:1:50;
[VgMesh, PoutMesh] = meshgrid(Vgrange, PoutRange);
[f,c] = contourf(Vgrange,PoutRange,F(VgMesh,PoutMesh),'ShowText','on');
xlabel('V_{g} (V)','FontSize',20,'FontName','Times New Roman');
ylabel('P_{out} (W)','FontSize',20,'FontName','Times New Roman');
title('Buck SC Power Loss','FontSize',24);
% ylim([0 80])

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',16)
colorbar
c.LevelList = [0 1 2 3 4 5 6 7 8 9];

hold on;
scatter(VgSim, PoutSim(locs),[],conditions(locs,2));


Vq = interp2(Vgrange,PoutRange,F(VgMesh,PoutMesh,,Yq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

stick = Ploss;
graph_values = 0;
return

mineff = .3;
etaSim(etaSim<mineff) = mineff;
etaSim(etaSim>1) = mineff;

locs = etaSim > mineff;

VgSim = conditions(locs,end);
F = scatteredInterpolant(VgSim, PoutSim(locs), PlossSim(locs), 'linear','none');


Vgrange = 16:.25:20;
PoutRange = 20:1:40;

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
stick = mean(mean(weighted_vals,'omitnan'),'omitnan');
stick = stick;

graph_values = F(X1,X2);
if isnan(stick)||isinf(stick)||stick<0
    stick = 100;
end

end

