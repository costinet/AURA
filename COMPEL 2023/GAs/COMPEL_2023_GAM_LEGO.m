clear

%%{
%% Set up optimization for 9 Variables
Nvars = 6;
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

LB = [1 1 1 1 10 0.3e6];
UB = [80 80 80 80 100 5e6];
sf = [1 1 1 1 0.1 1e-6];  %scale factor to ~ equalize magnitudes of variables

Aeq = [];
beq = [];
A = [];
b = [];

INTCON = [1 2 3 4];


%% Genetic Algorithm

options = optimoptions('gamultiobj','PlotFcn',{@gaplotpareto,@gaplotparetodistance,@gaplotrankhist,@gaplotspread},...
    'TolCon',1e-9,'Display','iter',...
    'MaxGenerations',75,'MaxTime',60*60*5, 'PopulationSize', 75);
options.InitialPopulation = [
    
   28.0000   48.0000   68.0000   52.0000    1.6905    0.3198
   28.0000   48.0000   49.0000   50.0000    3.3470    0.4593
   26.0000   14.0000   46.0000   13.0000    6.2559    1.0230
   28.0000   48.0000   68.0000   50.0000    4.6220    0.3858
   28.0000   48.0000   68.0000   51.0000    2.7808    0.3388
   28.0000   26.0000   48.0000   21.0000    6.2032    0.5894
   27.0000   15.0000   46.0000   14.0000    6.2559    1.0230
   28.0000   27.0000   61.0000   17.0000    5.4008    0.3405
   28.0000   48.0000   68.0000   51.0000    2.1076    0.3879
   28.0000   48.0000   68.0000   52.0000    1.7868    0.3225
   27.0000   46.0000   20.0000   13.0000    3.0520    0.5071
   28.0000   48.0000   68.0000   51.0000    2.4801    0.3642
   29.0000   48.0000   68.0000   50.0000    6.9382    0.3030
   28.0000   47.0000   20.0000   13.0000    3.0516    0.5071
   28.0000   16.0000   48.0000   17.0000    6.2076    0.6628
   28.0000   48.0000   68.0000   51.0000    3.7999    0.3872
   28.0000   48.0000   68.0000   52.0000    1.9408    0.3229
   27.0000   14.0000   46.0000   11.0000    6.4339    2.3469


    ];
% Add initial guess into initial population matrix to speed up process
[x,fval,exitflag,output,population,scores] = ...
    gamultiobj(@COMPEL_2023_LEGO_Sweep_GAM, Nvars, A, b, Aeq, beq, LB.*sf, UB.*sf, [], INTCON, options);
% save('gaResults.mat', 'x', 'fval', 'exitflag', 'output','population','scores')
% load('gaResults.mat');
%}


