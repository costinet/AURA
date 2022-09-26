function [stick,graph_values] = Texas_SC_Fib_Sweep(X,Coss_adj,Ron_adj)


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
% modelfile = 'AsyncBoost'; PLECsModel = 'Boost_Async';
% modelfile = 'MRBuck'; PLECsModel = 'MRBuck';
% modelfile = 'DSC4to1'; PLECsModel = 'HDSC';
% modelfile = 'DAB'; PLECsModel = 'DAB_oneCap';
% modelfile = 'DABfull'; PLECsModel = 'DAB_8Cap';
% modelfile = 'DSC4to1Diodes'; PLECsModel = 'HDSC_withDiodes';
%modelfile = 'SC_FIB_AURA_L.net';
%modelfile = 'Buck_Boost_Vout.net';
modelfile = 'SC_FIB_AURA.net';



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




%% Load for SC Fib


adjust = [ones(1,8), 1e-6, 100, 100] ;
X = X./adjust;

FETs = X(1:8);
fs = X(end-2);
dt1 = X(end-1);
dt2 = X(end);

% FET Selection
ron = [];
Coss = [];
FET_w = [];
FET_l = [];
for i = 1:length(FETs)
    [ron(i),Coss(i),FET_w(i),FET_l(i)]=Select_FET(FETs(i));
end

% Set switching interval
ON = 1;
OFF = 0;

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
    0     1     1     0     0     0     1     0     0     1     0     1     1     0     0     1
    0     1     1     0     1     0     1     0     0     1     0     1     1     0     0     1
    0     1     1     0     0     0     1     0     0     1     0     1     1     0     0     1
    0     1     1     1     0     1     1     0     0     1     0     1     1     1     0     1
    0     0     0     1     0     1     1     0     0     1     0     0     1     1     0     1
    1     0     0     1     0     1     1     0     0     1     0     0     1     1     0     1
    0     0     0     1     0     1     1     0     0     1     0     0     1     1     0     1
    ];

modSchemes(:,:,4) = [
    0     1     1     0     1     0     1     0     0     1     0     1     1     0     0     1
    0     1     1     0     0     0     1     0     0     1     0     1     1     0     0     0
    0     1     1     1     0     1     1     0     0     1     0     1     1     0     1     0
    0     0     0     1     0     1     1     0     0     1     0     0     1     0     0     0
    1     0     0     1     0     1     1     0     0     1     0     0     1     1     0     1
    0     0     0     1     0     1     1     0     0     1     0     0     0     1     0     1
    0     1     1     1     0     1     1     0     0     1     1     0     0     1     0     1
    0     1     1     0     0     0     1     0     0     1     0     0     0     0     0     1
    ];

modSchemes(:,:,5) = [
    0     1     1     1     0     1     1     0     0     1     0     1     1     0     1     0
    0     1     1     0     0     0     1     0     0     1     0     1     1     0     1     0
    0     1     1     0     1     0     1     0     0     1     0     1     1     0     1     0
    0     1     1     0     0     0     1     0     0     1     0     0     0     0     0     0
    0     1     1     1     0     1     1     0     0     1     1     0     0     1     0     1
    0     0     0     1     0     1     1     0     0     1     1     0     0     1     0     1
    1     0     0     1     0     1     1     0     0     1     1     0     0     1     0     1
    0     0     0     1     0     1     1     0     0     1     0     0     0     0     0     0
    ];

modSchemes(:,:,6) = [
    0     1     1     0     1     0     1     0     0     1     0     1     1     0     1     0
    0     1     1     0     0     0     1     0     0     1     0     0     0     0     1     0
    0     1     1     1     0     1     1     0     0     1     1     0     0     0     1     0
    0     0     0     1     0     1     1     0     0     1     1     0     0     0     0     0
    1     0     0     1     0     1     1     0     0     1     1     0     0     1     0     1
    0     0     0     1     0     1     1     0     0     1     1     0     0     0     0     0
    0     1     1     1     0     1     1     0     0     1     1     0     0     0     1     0
    0     1     1     0     0     0     1     0     0     1     0     0     0     0     1     0
    ];


swvec = [
    0     1     1     1     0     1     1     0     0     1     0     1     1     0     1     0
    0     1     1     0     0     0     1     0     0     1     0     1     1     0     1     0
    0     1     1     0     1     0     1     0     0     1     0     1     1     0     1     0
    0     1     1     0     0     0     1     0     0     1     0     0     0     0     0     0
    0     1     1     1     0     1     1     0     0     1     1     0     0     1     0     1
    0     0     0     1     0     1     1     0     0     1     1     0     0     1     0     1
    1     0     0     1     0     1     1     0     0     1     1     0     0     1     0     1
    0     0     0     1     0     1     1     0     0     1     0     0     0     0     0     0];


Vg = 12;


Vb1 = 4;
Vb2 = 4;

VgPOS = 1;
Vb1POS = 2;
Vb2POS = 3;



Rb = 5e-3;
RL = 20e-3;


Co = 2e-6; ESRo = 2e-3;
Cfly1 = 9.4e-6; ESR1 = 3e-3;  % 5V Cap
Cfly2 = 9.4e-6; ESR2 = 3e-3; % 10V Cap
Lc = 1.3e-6/2;
Ts = 1/fs;

u = [Vg Vb1 Vb2 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5]';


Numerical_Components = {
    'C1' Cfly1
    'C2' Cfly1
    'C3' Cfly2
    'C4' Cfly2
    'C5' Co
    'C6' Co
    'L1' Lc
    'L2' Lc
    'R1' Rb
    'R2' Rb
    'R3' RL
    'R4' RL
    'R5' ESR1
    'R6' ESR2
    'R7' ESR1
    'R8' ESR2
    'R9' ESRo
    'R10' ESRo
    'M1_C' Coss(1) + Coss_adj(1)
    'M2_C' Coss(2) + Coss_adj(2)
    'M3_C' Coss(3) + Coss_adj(3)
    'M4_C' Coss(2) + Coss_adj(2)
    'M5_C' Coss(1) + Coss_adj(1)
    'M6_C' Coss(3) + Coss_adj(3)
    'M7_C' Coss(4) + Coss_adj(4)
    'M8_C' Coss(5) + Coss_adj(5)
    'M9_C' Coss(5) + Coss_adj(5)
    'M10_C' Coss(4) + Coss_adj(4)
    'M11_C' Coss(6) + Coss_adj(6)
    'M12_C' Coss(7) + Coss_adj(7)
    'M13_C' Coss(8) + Coss_adj(8)
    'M14_C' Coss(7) + Coss_adj(7)
    'M15_C' Coss(6) + Coss_adj(6)
    'M16_C' Coss(8) + Coss_adj(8)
    'M1_R' ron(1) + Ron_adj(1)
    'M2_R' ron(2) + Ron_adj(2)
    'M3_R' ron(3) + Ron_adj(3)
    'M4_R' ron(2) + Ron_adj(2)
    'M5_R' ron(1) + Ron_adj(1)
    'M6_R' ron(3) + Ron_adj(3)
    'M7_R' ron(4) + Ron_adj(4)
    'M8_R' ron(5) + Ron_adj(5)
    'M9_R' ron(5) + Ron_adj(5)
    'M10_R' ron(4) + Ron_adj(4)
    'M11_R' ron(6) + Ron_adj(6)
    'M12_R' ron(7) + Ron_adj(7)
    'M13_R' ron(8) + Ron_adj(8)
    'M14_R' ron(7) + Ron_adj(7)
    'M15_R' ron(6) + Ron_adj(6)
    'M16_R' ron(8) + Ron_adj(8)
    };

ts = [0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25]./Ts;

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
    'M16_R'
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
    'M16'
    };



ON = 1;
OFF = 0;


% List all the resistances of the diodes or FETS when they are on or
% off
SW_OFF = ones(1,16).*10000000;

%SW_ON = [M1_R,M2_R,etc...]
SW_ON = [ron(1) ron(2) ron(3) ron(2) ron(1) ron(3) ron(4) ron(5) ron(5) ron(4) ron(6) ron(7) ron(8) ron(7) ron(6) ron(8)];
SW = [SW_OFF;SW_ON;SW_ON];


Diode_Forward_Voltage = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]'.*0;
u = [14.3 4 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
Order = [1 2 3 4 5 6 7 8];



etaSim = [];
PlossSim = [];
PoutSim = [];
conditions = [];
PossSim = [];
drange = linspace(0.05, .95, 12);
Vscloc1 = 32;
Vscloc2 = 33;
Illoc = 1;
Ib1loc = 2;
Ib2loc = 3;
VgPOS = 1;
i = 1;

u = u';

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


%% For the high input impedance circuit to test what Vg should be at 0W power
sim_Rb = SMPSim();
conv_Rb = sim_Rb.converter;
top_Rb = sim_Rb.topology;

% Set to run numerical parsesr

Numerical_Components(11,2) = {1e6};


%conv.ts = ts;
conv_Rb.u = u;
top_Rb.order = Order;
top_Rb.Element_Properties = Numerical_Components;
top_Rb.Switch_Resistors = Switch_Resistors;
top_Rb.Switch_Resistor_Values = SW;
top_Rb.Switch_Sequence = swvec;
top_Rb.Fwd_Voltage = Diode_Forward_Voltage;
top_Rb.Switch_Names = Switch_Names;
%
Numerical_Components(11,2) = {RL};




%% Analyze circuit
try
    for i = 4:size(modSchemes,3)
        
        swvec = modSchemes(:,:,i);
        top_Rb.Switch_Sequence = swvec;
        top_Rb.loadCircuit(modelfile,swvec,1);
        top.Switch_Sequence = swvec;
        top.loadCircuit(modelfile,swvec,1);
        
        for d = drange
            % dt = 0.0025;
            ds = [d-dt1, dt1, 1-d-dt2, dt2, d-dt1, dt1, 1-d-dt2, dt2];
            ts = ds/sum(ds)*Ts;
            
            conv_Rb.setSwitchingPattern(1:size(swvec,1), ts)
            % DMC_SteadyState(sim_Rb,conv_Rb,top_Rb);
            Xss = sim_Rb.steadyState;
            
            
            
            
            [ avgXs, avgYs ] = sim_Rb.ssAvgs(Xss);
            OLVin = avgYs(Vscloc1)+avgYs(Vscloc2);
            
            % MaxVin = OLVin+1;
            MaxVin = OLVin+0.5;
            if OLVin<14 || OLVin>21
                continue
            end
            Vinrange = OLVin:0.05:MaxVin;
            
            for Vin = Vinrange
                
               
                conv.setSwitchingPattern(1:size(swvec,1), ts)
                
                
                u(VgPOS) = Vin;
                sim.u = u;
                
                
                
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
                
                [valid]=DMC_SteadyState(sim,conv,top);
                if (~valid)
                    continue
                end
                Xss = sim.steadyState;
                
                [ avgXs, avgYs ] = sim.ssAvgs(Xss);
                
                
                Ib1 = (avgYs(47)+avgYs(23)*ESRo-Vb1)/Rb;
                Ib2 = (avgYs(48)+avgYs(24)*ESRo-Vb2)/Rb;
                I1 = -avgYs(17);
                
                eta = (Ib1*Vb1 + Ib2*Vb2)/(Vin*-I1 - Ib1^2*Rb - Ib2^2*Rb);
                Ploss = -(Ib1*Vb1 + Ib2*Vb2) + -(Vin*I1);
                Pout = (Ib1*Vb1 + Ib2*Vb2);
                
                
                if eta>1||eta<.3
                    continue
                end
                
                if Pout<0 || Ploss<0
                    continue
                end
                
                if Ploss>5
                    J = 445456465;
                end
                
                
                %for specific power
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

figure(104)

Vgrange = 14:.25:22;
PoutRange = 0:1:80;
[VgMesh, PoutMesh] = meshgrid(Vgrange, PoutRange);
[f,c] = contourf(Vgrange,PoutRange,F(VgMesh,PoutMesh),'ShowText','on');
xlabel('Vg');
ylabel('P_{out}');
% ylim([0 80])
colorbar
c.LevelList = [0 1 2 3 4 5 6 7 8 9];

hold on;
scatter(VgSim, PoutSim(locs),[],conditions(locs,2));

Vq = interp2(Vgrange,PoutRange,F(VgMesh,PoutMesh,,Yq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
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
stick = mean(mean(weighted_vals));
stick = stick;

graph_values = F(X1,X2);
if isnan(stick)||isinf(stick)||stick<0
    stick = 100;
end


end

