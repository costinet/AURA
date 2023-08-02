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

LB = [ 1   1  1  1  10   0.2e6];
UB = [ 80 80 64 10 200   5e6];
sf = [ 1   1  1  1   0.1 1e-6];  %scale factor to ~ equalize magnitudes of variables

Aeq = [];
beq = [];
A = [];
b = [];

INTCON = [1 2 3];


%% Genetic Algorithm

options = optimoptions('gamultiobj','PlotFcn',{@gaplotpareto,@gaplotparetodistance,@gaplotrankhist,@gaplotspread},...
    'TolCon',1e-9,'Display','iter',...
    'MaxGenerations',100,'MaxTime',60*60*1, 'PopulationSize', 100);
options.InitialPopulation = [
    
    7.0000   47.0000    9.0000  2  6.9939    2.0913
    8.0000   48.0000   35.0000  2  1.7978    0.5578
    8.0000   49.0000   13.0000  2  3.7423    0.6346
    8.0000   49.0000   34.0000  2  3.2075    0.7045
    6.0000   46.0000    9.0000  2  6.9939    2.0913
    8.0000   49.0000   38.0000  2  2.2510    0.4788
    8.0000   30.0000   38.0000  2  2.9308    0.5205
    8.0000   30.0000   38.0000  2  2.4605    0.4799
    8.0000   49.0000   35.0000  2  2.2776    0.6385
    8.0000   49.0000   35.0000  2  2.1475    0.5037
    8.0000   49.0000   14.0000  2  2.4514    0.7189
    8.0000   49.0000   34.0000  2  3.0783    0.6990
    8.0000   49.0000   35.0000  2  2.6541    0.6409
    8.0000   30.0000   38.0000  2  2.6345    0.4617
    8.0000   30.0000   38.0000  2  2.6853    0.5283
    8.0000   30.0000   35.0000  2  2.6689    0.5505
    7.0000   49.0000   29.0000  2  6.0766    1.1562
    7.0000   41.0000    9.0000  2  6.8805    1.1068

    ];
% Add initial guess into initial population matrix to speed up process
[x,fval,exitflag,output,population,scores] = ...
    gamultiobj(@COMPEL_2023_BuckHybridBuck_D_Sweep_GAM_CAP, Nvars, A, b, Aeq, beq, LB.*sf, UB.*sf, [], INTCON, options);
% save('gaResults.mat', 'x', 'fval', 'exitflag', 'output','population','scores')
% load('gaResults.mat');
%}



%{

%}






