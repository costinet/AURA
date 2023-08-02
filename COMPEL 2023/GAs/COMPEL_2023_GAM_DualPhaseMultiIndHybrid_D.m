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

LB = [ 1 1 1 10 0.3e6];
UB = [ 80 80 80 70 5e6];
sf = [ 1 1 1 0.1 1e-6];  %scale factor to ~ equalize magnitudes of variables

Aeq = [];
beq = [];
A = [];
b = [];

INTCON = [1 2 3];


%% Genetic Algorithm


options = optimoptions('gamultiobj','PlotFcn',{@gaplotpareto,@gaplotparetodistance,@gaplotrankhist,@gaplotspread},...
    'TolCon',1e-9,'Display','iter',...
    'MaxGenerations',50,'MaxTime',60*60*5, 'PopulationSize', 50);
options.InitialPopulation = [
    
    8.0000   29.0000   30.0000    2.4387    0.6853
   46.0000   49.0000   33.0000    5.9184    1.0002
   26.0000   46.0000   30.0000    6.3898    2.7960
   34.0000   46.0000   43.0000    6.8829    0.7448
   35.0000   49.0000   42.0000    6.1982    1.0887
    8.0000   29.0000   30.0000    3.5008    0.6845
   34.0000   47.0000   43.0000    6.9063    0.5035
   28.0000   49.0000   42.0000    6.8121    1.7095
   46.0000   50.0000   34.0000    6.1726    0.8825
   46.0000   50.0000   33.0000    5.9654    0.9603
    8.0000   30.0000   34.0000    3.0504    0.7118
    8.0000   29.0000   30.0000    3.0497    0.7217
    8.0000   49.0000   34.0000    4.2969    0.8512
    8.0000   30.0000   34.0000    3.1463    0.7323
    8.0000   29.0000   30.0000    2.8070    0.6966
   34.0000   49.0000   42.0000    6.8281    0.7768
    8.0000   29.0000   30.0000    2.7263    0.6949
   27.0000   47.0000   42.0000    6.7335    1.3192


    ];
% Add initial guess into initial population matrix to speed up process
[x,fval,exitflag,output,population,scores] = ...
    gamultiobj(@COMPEL_2023_DualPhaseMultiIndHybrid_D_Sweep_GAM, Nvars, A, b, Aeq, beq, LB.*sf, UB.*sf, [], INTCON, options);
% save('gaResults.mat', 'x', 'fval', 'exitflag', 'output','population','scores')
% load('gaResults.mat');
%}

