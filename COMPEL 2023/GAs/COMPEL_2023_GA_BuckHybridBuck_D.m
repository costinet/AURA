clear

%%{
%% Set up optimization for 9 Variables
Nvars = 5;
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

LB = [ 1 1 1 10 0.4e6];
UB = [ 80 80 80 40 3e6];
sf = [ 1 1 1 0.1 1e-6];  %scale factor to ~ equalize magnitudes of variables

Aeq = [];
beq = [];
A = [];
b = [];

INTCON = [1 2 3];

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
    'MaxGenerations',100,'MaxTime',60*60*1, 'MaxStallGenerations',50, 'PopulationSize', 30);
options.InitialPopulation = [
    
   27.0000   49.0000   29.0000    6.5562    1.7076
    8.0000   48.0000   38.0000    2.2031    0.6771
    8.0000   49.0000   34.0000    2.8232    0.7592
    9.0000   41.0000   34.0000    5.1239    0.6918
    7.0000   55.0000    9.0000    6.9971    2.1022
    8.0000   48.0000   35.0000    1.5753    0.4690
    7.0000   49.0000   34.0000    3.1276    0.7595
    6.0000   54.0000    9.0000    6.9971    2.1030
    7.0000   54.0000    9.0000    6.9971    2.1029
    7.0000   49.0000   34.0000    3.3948    0.6530
    8.0000   30.0000   38.0000    3.1598    0.5989
    8.0000   30.0000   38.0000    2.3871    0.4566
    8.0000   45.0000   34.0000    4.4334    0.7317
    8.0000   48.0000   30.0000    4.8117    0.9182
    8.0000   49.0000   34.0000    4.3195    0.7577
    6.0000   53.0000   10.0000    5.9836    1.5088
    8.0000   30.0000   38.0000    2.2464    0.4608
    6.0000   55.0000    9.0000    6.9971    2.1028


    ];
% Add initial guess into initial population matrix to speed up process
[x,fval,exitflag,output,population,scores] = ...
    ga(@COMPEL_2023_BuckHybridBuck_D_Sweep_GA, Nvars, A, b, Aeq, beq, LB.*sf, UB.*sf, [], INTCON, options);
% save('gaResults.mat', 'x', 'fval', 'exitflag', 'output','population','scores')
% load('gaResults.mat');
%}



%{

%}






