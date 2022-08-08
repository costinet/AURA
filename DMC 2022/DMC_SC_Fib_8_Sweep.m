function [stick,graph_values] = DMC_SC_Fib_8_Sweep(X,Coss_adj,Ron_adj)


niter = 0;
debug = 0;
debug2 = 0;

if ~debug
    w = warning ('off','all');
else
    w = warning ('on','all');
end

% sdir = mfilename('fullpath');
% sdir = sdir(1:find(sdir=='\',1,'last')-1);function [stick,graph_values] = DMC_SC_Fib_D_Sweep(X,Coss_adj,Ron_adj)


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
% modelfile = 'SC_FIB_AURA_L.net';
% modelfile = 'Buck_Boost_Vout.net';
% modelfile = 'SC_FIB_AURA_D.net';
% modelfile = 'SC_FIB_4_AURA.net';
modelfile = 'SC_FIB_8_AURA_D.net';

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


adjust = [ones(1,3), 1e-6, 100, 100] ;
X = X./adjust;

FETs = X(1:3);
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

ron(4) = ron(1);
ron(5) = ron(2);
ron(6) = ron(3);
Coss(4) = Coss(1);
Coss(5) = Coss(2);
Coss(6) = Coss(3);
Ron_adj(4) = Ron_adj(1);
Ron_adj(5) = Ron_adj(2);
Ron_adj(6) = Ron_adj(3);
Coss_adj(4) = Coss_adj(1);
Coss_adj(5) = Coss_adj(2);
Coss_adj(6) = Coss_adj(3);



% Set switching interval
ON = 1;
OFF = 0;
swvec = [
    OFF ON  ON  OFF ON  ON  OFF OFF OFF OFF OFF OFF   % Q1 ON
    OFF OFF OFF OFF OFF OFF OFF OFF OFF OFF OFF OFF  % Dead Q1 Q2
    ON  OFF OFF ON  OFF OFF OFF OFF OFF OFF OFF OFF % Q2 ON
    OFF OFF OFF OFF OFF OFF OFF OFF OFF OFF OFF OFF];  % Dead

Vg = 12;


Vb1 = 8;


VgPOS = 1;
Vb1POS = 2;
Vb2POS = 3;



Rb = 5e-3;
RL = 20e-3;


Co = 2e-6; ESRo = 2e-3;
Cfly1 = 9.4e-6; ESR1 = 3e-3;  % 5V Cap
Cfly2 = 9.4e-6; ESR2 = 3e-3; % 10V Cap
Lc = 2e-3/2;
Ts = 1/fs;

u = [Vg Vb1 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 ]';


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

Diode_Forward_Voltage = [0 0 0 0 0 0 1 1 1 1 1 1]'.*1.5;
u = [Vg Vb1 0 0 0 0 0 0 1.5 1.5 1.5 1.5 1.5 1.5]';
Order = [1 2 3 4];



etaSim = [];
PlossSim = [];
PoutSim = [];
conditions = [];
PossSim = [];
drange = linspace(0.15, .65, 8);
Vscloc1 = 48;
Vscloc2 = 49;
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
    for i = 1
        
        
        top_Rb.Switch_Sequence = swvec;
        top_Rb.loadCircuit(modelfile,swvec,1);
        top.Switch_Sequence = swvec;
        top.loadCircuit(modelfile,swvec,1);
        
        for d = drange
            % dt = 0.0025;
            ds = [d-dt1, dt1, 1-d-dt2, dt2];
            ts = ds/sum(ds)*Ts;
            
            conv_Rb.setSwitchingPattern(1:size(swvec,1), ts)
            % DMC_SteadyState(sim_Rb,conv_Rb,top_Rb);
            Xss = sim_Rb.steadyState;
            
            if(debug)
                    sim_Rb.plotAllStates(1);
                end
            
            
            [ avgXs, avgYs ] = sim_Rb.ssAvgs(Xss);
            %OLVin = avgYs(25)+avgYs(26)+avgYs(27);
            OLVin = ((1-d)*16)+8;
            % MaxVin = OLVin+1;
            MaxVin = OLVin+0.5;
            %if OLVin<14 || OLVin>21
            %    continue
            %end
            Vinrange = OLVin:0.1:MaxVin;
            
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
                
                
                Ib1 = (avgYs(24)+avgYs(7)*ESRo-Vb1)/Rb;
               
                I1 = -avgYs(4);
                
                eta = (Ib1*Vb1 )/(Vin*-I1 - Ib1^2*Rb );
                Ploss = -(Ib1*Vb1) + -(Vin*I1);
                Pout = (Ib1*Vb1 );
                
                
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

figure(106)

Vgrange = 14:.25:22;
PoutRange = 0:1:50;
[VgMesh, PoutMesh] = meshgrid(Vgrange, PoutRange);
[f,c] = contourf(Vgrange,PoutRange,F(VgMesh,PoutMesh),'ShowText','on');
xlabel('V_{g} (V)','FontSize',20,'FontName','Times New Roman');
ylabel('P_{out} (W)','FontSize',20,'FontName','Times New Roman');
title('SC Fib 8V Power Loss','FontSize',24);
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

