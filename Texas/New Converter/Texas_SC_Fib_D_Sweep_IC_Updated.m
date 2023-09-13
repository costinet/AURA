function [stick,graph_values] = Texas_SC_Fib_D_Sweep_IC_Updated(X)


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

modelfile = 'SC_FIB_TI_AURA_D.net';


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
Ron_adj = zeros(13,1);
Coss_adj = zeros(13,1);
adjust = [ones(1,13).*20, 1e-6] ;
X = X./adjust;

FETs = X(1:13);
fs = X(end);


% FET Selection
ron = [];
Coss = [];
FET_w = [];
FET_l = [];
tot_area = [];

W_selection = FETs;


for select_FET = 1:length(W_selection)
            a =  0.001447;
            b = -1;
            ron(select_FET) = a*W_selection(select_FET)^b;
            p1 = 2.487e-9;
            p2 = -6.393e-13;
            Coss(select_FET) = p1*W_selection(select_FET)+p2;
            
            L3=900*10^-9;
            
            tot_area = tot_area+W_selection(select_FET)*L3;
end

% Set switching interval
ON = 1;
OFF = 0;

modSchemes(:,:,1) = [
    0 1 1 0 1 1 1 1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 1 1 0 1 1 1 1 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 1 1 0 1 1 1 1 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 % needed to fix size problem
];


modSchemes(:,:,2) = [
    0 1 1 0 1 1 1 1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 1 1 1 0 0 1 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 1 1 1 0 0 1 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0% needed to fix size problem 
];


modSchemes(:,:,3) = [
    0 1 1 0 1 0 0 1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 1 1 1 0 0 1 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0
    1 0 0 0 1 1 1 1 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0
]; 

modSchemes(:,:,4) = [
    1 0 0 0 1 1 1 1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 1 1 0 1 0 0 1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 1 1 1 0 0 1 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0
];

modSchemes(:,:,5) = [
    1 0 0 1 0 0 1 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 1 1 0 1 0 0 1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    1 0 0 0 1 1 1 1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
];

modSchemes(:,:,6) = [
    1 0 0 0 1 1 1 1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 1 1 0 1 0 0 1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    1 0 0 1 0 0 1 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0
];


swvec = [
    1 0 0 0 1 1 1 1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 1 1 0 1 0 0 1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    1 0 0 1 0 0 1 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0
       ];


Vg = 12;


Vb1 = 4;
Vb2 = 4;

VgPOS = 1;
Vb1POS = 2;
Vb2POS = 3;



Rb = 5e-3;
RL = 320e-3/2;


Co = 20e-6; ESRo = 2e-3;
Cfly1 = 9.4e-6; ESR1 = 3e-3;  % 5V Cap
Cfly2 = 9.4e-6; ESR2 = 3e-3; % 10V Cap
Lc = 1.3e-6/2;
fs = 1256200;
Ts = 1/fs;

% This is from expierimental testing that occured on in early 2023 that
% tried to determine what the acutal values are of the circuit. So far the
% Ron values really hurt the performance of the converter but are accurate
% to the acutal values used in the performance


%Coss = [4.5e-10 6e-10 6e-10 6e-10 6e-10 ];
%ron  = [0.1 0.1 0.1 0.1 0.1]; 

%u = [Vg Vb1 Vb2 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5]';


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
    'R11' ESRo
    'R12' ESRo
    'M1_C' Coss(1) + Coss_adj(1)
    'M2_C' Coss(2) + Coss_adj(2)
    'M3_C' Coss(3) + Coss_adj(3)
    'M4_C' Coss(4) + Coss_adj(4)
    'M5_C' Coss(5) + Coss_adj(4)
    'M6_C' Coss(6) + Coss_adj(1)
    'M7_C' Coss(7) + Coss_adj(2)
    'M8_C' Coss(8) + Coss_adj(3)
    'M9_C' Coss(9) + Coss_adj(4)
    'M10_C' Coss(10) + Coss_adj(5)
    'M11_C' Coss(11) + Coss_adj(11)
    'M12_C' Coss(12) + Coss_adj(12)
    'M13_C' Coss(13) + Coss_adj(13)
    'D1_C' Coss(1) + Coss_adj(1)
    'D2_C' Coss(2) + Coss_adj(2)
    'D3_C' Coss(3) + Coss_adj(3)
    'D4_C' Coss(4) + Coss_adj(4)
    'D5_C' Coss(5) + Coss_adj(5)
    'D6_C' Coss(6) + Coss_adj(6)
    'D7_C' Coss(7) + Coss_adj(7)
    'D8_C' Coss(8) + Coss_adj(8)
    'D9_C' Coss(9) + Coss_adj(9)
    'D10_C' Coss(10) + Coss_adj(10)
    'D11_C' Coss(11) + Coss_adj(11)
    'D12_C' Coss(12) + Coss_adj(12)
    'D13_C' Coss(13) + Coss_adj(13)
    'M1_R' ron(1) + Ron_adj(1)
    'M2_R' ron(2) + Ron_adj(2)
    'M3_R' ron(3) + Ron_adj(3)
    'M4_R' ron(4) + Ron_adj(4)
    'M5_R' ron(5) + Ron_adj(5)
    'M6_R' ron(6) + Ron_adj(6)
    'M7_R' ron(7) + Ron_adj(7)
    'M8_R' ron(8) + Ron_adj(8)
    'M9_R' ron(9) + Ron_adj(9)
    'M10_R' ron(10) + Ron_adj(10)
    'M11_R' ron(11) + Ron_adj(11)
    'M12_R' ron(12) + Ron_adj(12)
    'M13_R' ron(13) + Ron_adj(13)
    'D1_R' ron(1) + Ron_adj(1)
    'D2_R' ron(2) + Ron_adj(2)
    'D3_R' ron(3) + Ron_adj(3)
    'D4_R' ron(4) + Ron_adj(4)
    'D5_R' ron(5) + Ron_adj(5)
    'D6_R' ron(6) + Ron_adj(6)
    'D7_R' ron(7) + Ron_adj(7)
    'D8_R' ron(8) + Ron_adj(8)
    'D9_R' ron(9) + Ron_adj(9)
    'D10_R' ron(10) + Ron_adj(10)
    'D11_R' ron(11) + Ron_adj(11)
    'D12_R' ron(12) + Ron_adj(12)
    'D13_R' ron(13) + Ron_adj(13)

    };

ts = [0.25 0.25 ]./Ts;

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
    };



ON = 1;
OFF = 0;


% List all the resistances of the diodes or FETS when they are on or
% off
SW_OFF = ones(1,26).*10000000;

%SW_ON = [M1_R,M2_R,etc...]
SW_ON = [
    ron(1)+Ron_adj(1) 
    ron(2)+Ron_adj(2) 
    ron(3)+Ron_adj(3) 
    ron(4)+Ron_adj(4)
    ron(5)+Ron_adj(5) 
    ron(6)+Ron_adj(6)
    ron(7)+Ron_adj(7)
    ron(8)+Ron_adj(8)
    ron(9)+Ron_adj(9)
    ron(10)+Ron_adj(10)
    ron(11)+Ron_adj(11)
    ron(12)+Ron_adj(12) 
    ron(13)+Ron_adj(13)
    ron(1)+Ron_adj(1) 
    ron(2)+Ron_adj(2) 
    ron(3)+Ron_adj(3) 
    ron(4)+Ron_adj(4)
    ron(5)+Ron_adj(5) 
    ron(6)+Ron_adj(6)
    ron(7)+Ron_adj(7)
    ron(8)+Ron_adj(8)
    ron(9)+Ron_adj(9)
    ron(10)+Ron_adj(10)
    ron(11)+Ron_adj(11)
    ron(12)+Ron_adj(12) 
    ron(13)+Ron_adj(13)
    ];
SW = [SW_OFF;SW_ON';SW_ON'];


Diode_Forward_Voltage = [ 0 0 0 0 0 0 0 0 0 0 0 0 0   1 1 1 1 1 1 1 1 1 1 1 1 1]'.*0.61;
% u = [14.3 4 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5]';
u = [14.3 4 0 0 0 0 0 0 0 0 0 0 0 0  0  0.61 0.61 0.61 0.61 0.61 0.61 0.61 0.61 0.61 0.61 0.61 0.61 0.61 ]';
Order = [1 2 3];



etaSim = [];
PlossSim = [];
PoutSim = [];
conditions = [];
PossSim = [];
drange = linspace(0.15, .85, 6);
Vscloc1 = 41;
Vscloc2 = 42;
Vscloc3 = 45;

%Vscloc1 = 50;
%Vscloc2 = 51;
%Vscloc1 = 56;
%0Vscloc2 = 57;
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
    for i = 1:6
        if i==3
            J = 456456465;
        end
        swvec =modSchemes(:,:,i);
        top_Rb.Switch_Sequence = swvec;
        top_Rb.loadCircuit(modelfile,swvec,1);
        top.Switch_Sequence = swvec;
        top.loadCircuit(modelfile,swvec,1);

        for d = drange
            % dt = 0.0025;
            ds = [d 1-d-0.05 0.05];
            ts = ds/sum(ds)*Ts;
            
            conv_Rb.setSwitchingPattern(1:size(swvec,1), ts)
            % DMC_SteadyState(sim_Rb,conv_Rb,top_Rb);
            Xss = sim_Rb.steadyState;
            

            [ avgXs, avgYs ] = sim_Rb.ssAvgs(Xss);
            OLVin = avgYs(Vscloc1)+avgYs(Vscloc2)+avgYs(Vscloc3);
            
            % MaxVin = OLVin+1;
                MaxVin = OLVin+5;

            if OLVin<2 || OLVin>22
                continue
            end
            Vinrange = OLVin:0.4:MaxVin;
            
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
                
                %[valid]=DMC_SteadyState(sim,conv,top);

%                 if (~valid)
%                     continue
%                 end
                Xss = sim.steadyState;
                
                [ avgXs, avgYs ] = sim.ssAvgs(Xss);
                
                [ xs, t, ys] = sim.SS_WF_Reconstruct;
                



                Ib1 = (avgYs(62)+avgYs(30)*ESRo-Vb1)/Rb;
               
                I1 = -avgYs(27);
                
                eta = (Ib1*Vb1 )/(Vin*-I1 - Ib1^2*Rb);
                Ploss = -(Ib1*Vb1) + -(Vin*I1);
                Pout = (Ib1*Vb1);
                
                
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
                
                    if i == 3 && d < 0.45 && d>0.4 
                        J = 32423423;
                    end


                if Pout > 40
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

Vgrange = 5:.25:22;
PoutRange = 0:1:50;
[VgMesh, PoutMesh] = meshgrid(Vgrange, PoutRange);
[f,c] = contourf(Vgrange,PoutRange,F(VgMesh,PoutMesh),'ShowText','on');
xlabel('V_{g} (V)','FontSize',20,'FontName','Times New Roman');
ylabel('P_{out} (W)','FontSize',20,'FontName','Times New Roman');
title('SC Fib Power Loss','FontSize',24);
% ylim([0 80])

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',16)
colorbar
c.LevelList = [0 1 2 3 4 5 6 7 8 9];

hold on;
scatter(VgSim, PoutSim(locs),[],conditions(locs,2));

figure(105)
scatter(VgSim, PoutSim(locs),[],conditions(locs,2));

Vq = interp2(Vgrange,PoutRange,F(VgMesh,PoutMesh),Xq6,Yq6)
Vq = interp2(Vgrange,PoutRange,F(VgMesh,PoutMesh),Xq,Yq)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
mineff = .3;
etaSim(etaSim<mineff) = mineff;
etaSim(etaSim>1) = mineff;

locs = etaSim > mineff;

VgSim = conditions(locs,end);
F = scatteredInterpolant(VgSim, PoutSim(locs), PlossSim(locs), 'linear','none');


Vgrange = 7:1:20;
PoutRange = 20:1:30;

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

