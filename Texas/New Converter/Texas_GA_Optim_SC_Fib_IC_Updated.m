

clear

%%{
%% Set up optimization for 9 Variables
Nvars = 14;
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


LB = [ 0.005 0.005 0.005 0.005 0.005 0.005 0.005 0.005 0.005 0.005 0.005 0.005 0.005 0.9e6];
UB = [ 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5  1.5e6];
sf = [ 20 20 20 20 20 20 20 20 20 20 20 20 20  1e-6];  %scale factor to ~ equalize magnitudes of variables


L23=900*10^-9;

Aeq = [1 1 1 1 1 1 1 1 1 1 1 1 1 0].*L23./sf;
beq = 3e-6;


%Aeq = [];
%beq = [];
A = [];
b = [];
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

options = optimoptions('ga','PlotFcn',@gaplotbestf,'Display','iter',...
    'MaxGenerations',500,'MaxTime',60*60*9, 'MaxStallGenerations',50,...
    'PopulationSize', 100,'ConstraintTolerance',1e-10);
options.InitialPopulation = [

0.3803    0.2511    0.4304    0.2612    0.3086    0.2328    0.3234    0.2331    0.1130    0.2691    0.1690    0.2551    0.1063   1.1698e6

].*sf;
% Add initial guess into initial population matrix to speed up process
[x,fval,exitflag,output,population,scores] = ...
    ga(@Texas_SC_Fib_D_Sweep_IC_Updated, Nvars, A, b, Aeq, beq, LB.*sf, UB.*sf, [], INTCON, options);
% save('gaResults.mat', 'x', 'fval', 'exitflag', 'output','population','scores')
% load('gaResults.mat');
%}

x
fval
exitflag
output
population
scores




%{
First run with unlimited area contraint M1-13 plus f 30-Aug-2023 12:15:20
0.3941    0.3652    0.3978    0.3631    0.3694    0.3269    0.3840    0.3387    0.1661    0.3324    0.2045    0.3399    0.1601  1.2520e6

Next run without constraint by accident
0.4408    0.3652    0.3978    0.3631    0.3694    0.3269    0.3840   0.3382    0.1661    0.3324    0.2045    0.3399    0.1601   1.7716e6

Run with 3e-6 area constraint but I believe tolerances are a bit off  30-Aug-2023 18:45:41
0.4458    0.3528    0.4025    0.3684    0.3619    0.3195    0.4146   0.3332    0.1584    0.3884    0.1958    0.3768    0.1513   1.1612e6

Run with 3e-6 area constraint where it was acutally enforced 30-Aug-2023 20:14:32
0.3892    0.2547    0.4369    0.2671    0.2935    0.2390    0.2997    0.2350    0.1185    0.2753    0.1545    0.2600    0.1101
0.389174271969520   0.254651448393115   0.436866270499415   0.267097607771296   0.293521151300079   0.239011891932148   0.299727595809002   0.235034512607996   0.118508195326722   0.275300430127527   0.154455153335458   0.259981748398014   0.110069722529707

% Run overnight with beq = 3e-6  31-Aug-2023 14:40:12
0.3803    0.2511    0.4304    0.2612    0.3086    0.2328    0.3234    0.2331    0.1130    0.2691    0.1690    0.2551    0.1063   1.1698e6





%}

