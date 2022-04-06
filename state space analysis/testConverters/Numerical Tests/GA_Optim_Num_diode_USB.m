

clear

%%{
%% Set up optimization for 9 Variables
Nvars = 11;
% Variables in order:
%{
 FET Selection:
M1 and M5
M2 and M4
M3 and M6
M7 and M10
M8 and M9
M11 and M15
M12 and M14
M13 and M16

FET Width:
M1 and M5
M2 and M4
M3 and M6
M7 and M10
M8 and M9
M11 and M15
M12 and M14
M13 and M16

frequency

%}

LB = [0.005 0.005 0.005 0.005 0.005 0.005 0.005 0.005 1e6 0.001 0.001];
UB = [0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 2e6 0.0075 0.0075];
sf = [20 20 20 20 20 20 20 20 1e-6 100 100];  %scale factor to ~ equalize magnitudes of variables


L23=900*10^-9;

Aeq = [];
beq = [];
A = [2 2 2 2 2 2 2 2 0 0 0].*L23./sf;
b = 1.0001e-6;

INTCON = [];

%% Particle Swarm
% options = optimoptions(@particleswarm,'PlotFcn','pswplotbestf');
% [x,fval,exitflag,output] = particleswarm(@DABLossMimize, Nvars, LB.*sf, UB.*sf, options);
% save('psResults.mat', 'x', 'fval', 'exitflag', 'output')
% load('psResults.mat');

%% Genetic Algorithm

% Constraints
% These aren't necessary, but can speed up the process by eliminating
% regions of the design space that we can conclude outright won't meet our
% goal.  In this case, neglect anything worse than the 97.5% design we
% already have, based on the ideal model
%A = [Ig^2*2, Iout^2*2/1000, 0,0, Ig^2/Q1/1e6]; %require ideal conduction loss less than 2.5% of Pout
%b = Pout*.025;

options = optimoptions('ga','PlotFcn',@gaplotbestf,'TolCon',1e-9,'Display','iter',...
    'MaxGenerations',100,'MaxTime',60*60*4, 'MaxStallGenerations',50, 'PopulationSize', 20);
options.InitialPopulation = [

    0.1355    0.0843    0.2507    0.3038    0.0274    0.0764    0.0461    0.1828    1.3327e6    0.004395  0.001267
    0.1389    0.0793    0.2542    0.3055    0.0271    0.0763    0.0424    0.1863    1.2951e6    0.004434  0.001225
    0.1042    0.0595    0.1906    0.2291    0.0203    0.0572    0.0318    0.1397    1.2951e6    0.004434  0.001225
    0.0694    0.0396    0.1271    0.1527    0.0135    0.0382    0.0212    0.0931    1.2951e6    0.004434  0.001225
    
    ].*sf;
% Add initial guess into initial population matrix to speed up process
[x,fval,exitflag,output,population,scores] = ...
    ga(@AURA_Eff_Sweep_IC_reduce_USB, Nvars, A, b, Aeq, beq, LB.*sf, UB.*sf, [], INTCON, options);
% save('gaResults.mat', 'x', 'fval', 'exitflag', 'output','population','scores')
% load('gaResults.mat');
%}



%% Saved Results
%{
0.1273    0.0608    0.2518    0.3007    0.0050    0.1087    0.0474    0.2074    1.5244e6    0.0045    0.0010
    0.1268    0.0618    0.2513    0.2994    0.0050    0.1086    0.0505    0.2078    1.5250e6    0.0045    0.0010
    0.1285    0.0599    0.2469    0.2988    0.0050    0.0928    0.0477    0.2311    1.5177e6    0.0045    0.0010
    0.1285    0.0641    0.2768    0.3030    0.0050    0.0871    0.0524    0.1885    1.5480e6    0.0045    0.0010
    0.1295    0.0674    0.2780    0.2988    0.0050    0.0842    0.0588    0.1893    1.5902e6    0.0045    0.0010
    0.1097    0.0686    0.2877    0.2794    0.0052    0.0880    0.0590    0.2125    1.5662e6    0.0045    0.0010
    0.1217    0.0669    0.2847    0.2795    0.0051    0.0881    0.0572    0.2082    1.5709e6    0.0045    0.0010
    0.1741    0.0758    0.2714    0.3290    0.0415    0.0924    0.0530    0.2747    1.2015e6    0.0049    0.0010
    0.1859    0.0751    0.2862    0.3284    0.0425    0.0920    0.0525    0.2743    1.2015e6    0.0049    0.0010
    0.1885    0.0816    0.3029    0.3431    0.0353    0.0932    0.0501    0.2841    1.2128e6    0.0059    0.0013

%}
