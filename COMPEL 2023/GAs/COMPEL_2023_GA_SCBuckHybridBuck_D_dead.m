clear

%%{
%% Set up optimization for 9 Variables
Nvars = 8;
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

LB = [1 1 1 1 1 1 10 0.4e6];
UB = [80 80 80 80 80 80 50 3.5e6];
sf = [1 1 1 1 1 1 0.1 1e-6];  %scale factor to ~ equalize magnitudes of variables

Aeq = [];
beq = [];
A = [];
b = [];

INTCON = [1 2 3 4 5 6];

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
    'MaxGenerations',1000,'MaxTime',60*60*6, 'MaxStallGenerations',50, 'PopulationSize', 20);
options.InitialPopulation = [
    
       17 8 1 1 1 6 10 2.108e6 
       6 7 65 13 3 64 17.4 2.0950e6

    ].*sf;
% Add initial guess into initial population matrix to speed up process
[x,fval,exitflag,output,population,scores] = ...
    ga(@COMPEL_2023_SCBuckHybridBuck_D_Sweep_dead_GA, Nvars, A, b, Aeq, beq, LB.*sf, UB.*sf, [], INTCON, options);
% save('gaResults.mat', 'x', 'fval', 'exitflag', 'output','population','scores')
% load('gaResults.mat');
%}



%{

%}






