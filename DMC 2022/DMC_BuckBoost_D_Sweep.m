function [stick,graph_values] = DMC_BuckBoost_Sweep(X,Coss_adj,Ron_adj)


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
modelfile = 'Buck_Boost_Vout_D.net';
% modelfile = 'Buck_Boost_Vout.net';
%modelfile = 'SC_FIB_AURA.net';



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




%% Load for 3 Level Buck

%%{
adjust = [ones(1,6), 1, 1e-6, 100, 100] ;
X = X./adjust;

FETs = X(1:6);
Inductor_Chosen = X(end-3);
fs = X(end-2);
dt1 = X(end-1);
dt2 = X(end);


% Inductor Selection
Selected_Inductor = Select_Inductor(Inductor_Chosen);

R1_L  = Selected_Inductor(1);
R2_L = Selected_Inductor(2);
C_L = Selected_Inductor(3);
k1 = Selected_Inductor(4);
k2  = Selected_Inductor(5);
k3 = Selected_Inductor(6);
k4 = Selected_Inductor(7);
k5 = Selected_Inductor(8);
L_area = Selected_Inductor(9);

C_L = C_L*1e-12;
Rvar1 = k1*sqrt(fs);
Rvar2 = k2*sqrt(fs);
Lvar = (k3-k4*log(k5*fs))*1e-6;

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
swvec = [
    ON  OFF ON  ON  OFF OFF OFF OFF OFF OFF OFF OFF   % Q1 ON
    OFF OFF ON  ON  OFF OFF OFF OFF OFF OFF OFF OFF  % Dead Q1 Q2
    OFF ON  ON  ON  OFF OFF OFF OFF OFF OFF OFF OFF % Q2 ON
    OFF OFF ON  ON  OFF OFF OFF OFF OFF OFF OFF OFF];  % Dead


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

Cfly = 23.5e-6;
CflyESR = 2.2e-3;

Vb1 = Vbat1;
Vb2 = Vbat2;
Rb = Rbatt;



Ts = 1/fs;
dt = 5e-9;

Numerical_Components = {
    'C1' Cfly
    'C2' Cout
    'C3' Cout
    'L1' Lvar
    'M1_C' Coss(1)+Coss_adj(1)
    'M2_C' Coss(2)+Coss_adj(2)
    'M3_C' Coss(3)+Coss_adj(3)
    'M4_C' Coss(4)+Coss_adj(4)
    'M5_C' Coss(5)+Coss_adj(5)
    'M6_C' Coss(6)+Coss_adj(6)
    'D1_C' Coss(1)+Coss_adj(1)
    'D2_C' Coss(2)+Coss_adj(2)
    'D3_C' Coss(3)+Coss_adj(3)
    'D4_C' Coss(4)+Coss_adj(4)
    'D5_C' Coss(5)+Coss_adj(5)
    'D6_C' Coss(6)+Coss_adj(6)
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
    'M1_R' ron(1)+Ron_adj(1)
    'M2_R' ron(2)+Ron_adj(2)
    'M3_R' ron(3)+Ron_adj(3)
    'M4_R' ron(4)+Ron_adj(4)
    'M5_R' ron(5)+Ron_adj(5)
    'M6_R' ron(6)+Ron_adj(6)
    'D1_R' ron(1)+Ron_adj(1)
    'D2_R' ron(2)+Ron_adj(2)
    'D3_R' ron(3)+Ron_adj(3)
    'D4_R' ron(4)+Ron_adj(4)
    'D5_R' ron(5)+Ron_adj(5)
    'D6_R' ron(6)+Ron_adj(6)
    };

ts_Baxter = [0.70 0.001 0.3 0.001];
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
    'D1_R'
    'D2_R'
    'D3_R'
    'D4_R'
    'D5_R'
    'D6_R'
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
    'D1'
    'D2'
    'D3'
    'D4'
    'D5'
    'D6'
    };



ON = 1;
OFF = 0;

% List all the resistances of the diodes or FETS when they are on or
% off
SW_OFF = ones(1,12).*10000000;

%SW_ON = [M1_R,M2_R,etc...]
SW_ON = [ron(1)+Ron_adj(1) ron(2)+Ron_adj(2) ron(3)+Ron_adj(3) ron(4)+Ron_adj(4) ron(5)+Ron_adj(5) ron(6)+Ron_adj(6) ron(1)+Ron_adj(1) ron(2)+Ron_adj(2) ron(3)+Ron_adj(3) ron(4)+Ron_adj(4) ron(5)+Ron_adj(5) ron(6)+Ron_adj(6)];
SW = [SW_OFF;SW_ON;SW_ON];

Diode_Forward_Voltage = [1 1 1 1 1 1 1 1 1 1 1 1]'.*1;
u = [Vg Vbat1 Vbat2 0 0 0 0 0 0 1 1 1 1 1 1]';
Order = [1 2 3 4];
%}
etaSim = [];
PlossSim = [];
PoutSim = [];
conditions = [];
PossSim = [];
drange = linspace(0.39, .6, 12);
Vscloc1 = 15;
Vscloc2 = 18;
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




%% Analyze circuit
try
    for d = drange
        % dt = 0.0025;
        ds = [d-dt1, dt1, 1-d-dt2, dt2];
        ts = ds/sum(ds)*Ts;
        
        OLVin = 8/d;
        % MaxVin = OLVin+1;
        MaxVin = OLVin+2.5;
        if OLVin<14 || OLVin>21
            continue
        end
        Vinrange = OLVin:.05:MaxVin;
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
            
            
            Ib1 = (avgYs(27)+avgYs(10)*Rcap-Vbat1)/Rbatt;
            Ib2 = (avgYs(25)+avgYs(8)*Rcap-Vbat2)/Rbatt;
            I1 = -avgYs(1);
            
            eta = (Ib1*Vb1 + Ib2*Vb2)/(Vin*-I1 - Ib1^2*Rb - Ib2^2*Rb);
            Ploss = -(Ib1*Vb1 + Ib2*Vb2) + -(Vin*I1);
            Pout = (Ib1*Vb1 + Ib2*Vb2);
            
            
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

