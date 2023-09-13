%{
 __      ___ _   _      ___              _ _   _           
 \ \    / (_) |_| |_   |   \ ___ __ _ __| | |_(_)_ __  ___ 
  \ \/\/ /| |  _| ' \  | |) / -_) _` / _` |  _| | '  \/ -_)
   \_/\_/ |_|\__|_||_| |___/\___\__,_\__,_|\__|_|_|_|_\___|
%}                                                         
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

A = [1 1 1 1 1 1 1 1 1 1 1 1 1 0].*L23./sf;
b = 2.01e-6;


Aeq = [];
beq = [];
%A = [];
%b = [];
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
 
  % 0.3803    0.2511    0.4304    0.2612    0.3086    0.2328    0.3234    0.2331    0.1130    0.2691    0.1690    0.2551    0.1063    1.1698e6
  % 0.3849    0.2431    0.4338    0.2636    0.3094    0.2368    0.3180    0.2343    0.1152    0.2713    0.1741    0.2535    0.1090    1.2132e6
  % 0.3959    0.2388    0.4233    0.2686    0.3369    0.2309    0.3159    0.2407    0.1200    0.2097    0.1613    0.2536    0.1148    1.1708e6
  % 0.4131    0.2422    0.4340    0.2280    0.3222    0.2342    0.3153    0.2310    0.1209    0.2740    0.1630    0.2485    0.1101    1.3131e6
  %  0.4171    0.4975    0.4624    0.2101    0.3651    0.2412    0.3059    0.2124    0.0704    0.2093    0.0652    0.2159    0.0580    1.3026e6
  %  0.2781    0.3317    0.3083    0.1401    0.2434    0.1608    0.2039    0.1416    0.0469    0.1395    0.0435    0.1439    0.0387    1.3026e6
    0.2971    0.2414    0.3070    0.1384    0.2967    0.2368    0.1602    0.1087    0.0423    0.1316    0.0611    0.1366    0.0560    1.4088e6

].*sf;
% Add initial guess into initial population matrix to speed up process
[x,fval,exitflag,output,population,scores] = ...
    ga(@Texas_SC_Fib_D_Sweep_IC_Updated_dead, Nvars, A, b, Aeq, beq, LB.*sf, UB.*sf, [], INTCON, options);
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

% Run during day same as above with b<3.001e-6 31-Aug-2023 19:51:54 with
slight improvement from above
0.3849    0.2431    0.4338    0.2636    0.3094    0.2368    0.3180    0.2343    0.1152    0.2713    0.1741    0.2535    0.1090    1.2132e6

% Run overnight with dt = 0.0005 and updated 12-16 V modulation scheme for
improved eff and no overlap  01-Sep-2023 13:27:34
0.3959    0.2388    0.4233    0.2686    0.3369    0.2309    0.3159    0.2407    0.1200    0.2097    0.1613    0.2536    0.1148    1.1708e6

% Run overnight with same as above 06-Sep-2023 10:51:52
0.4131    0.2422    0.4340    0.2280    0.3222    0.2342    0.3153    0.2310    0.1209    0.2740    0.1630    0.2485    0.1101    1.3131e6

% Run overnight with reduced fied of fval focused on 20 V and 30ish Watts 07-Sep-2023 11:32:44
0.4171    0.4975    0.4624    0.2101    0.3651    0.2412    0.3059    0.2124    0.0704    0.2093    0.0652    0.2159    0.0580    1.3026e6

% Short run to confirm above results 08-Sep-2023 19:56:23
% Need to reduce total area constraint to accomidate hitting the boundary on M2
0.4401    0.4950    0.4599    0.2326    0.3638    0.2755    0.3062    0.1974    0.0398    0.2109    0.0634    0.2151    0.0424    1.3026e6

% Run at  09-Sep-2023 13:54:04
0.2971    0.2414    0.3070    0.1384    0.2967    0.2368    0.1602    0.1087    0.0423    0.1316    0.0611    0.1366    0.0560    1.4088e6



%}

