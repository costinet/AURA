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

LB = [  1  1  1 10    0.3e6];
UB = [ 80 80 80 70    5e6];
sf = [  1  1  1  0.1  1e-6];  %scale factor to ~ equalize magnitudes of variables

Aeq = [];
beq = [];
A = [];
b = [];

INTCON = [1 2 3];


%% Genetic Algorithm

options = optimoptions('gamultiobj','PlotFcn',{@gaplotpareto,@gaplotparetodistance,@gaplotrankhist,@gaplotspread},...
    'TolCon',1e-9,'Display','iter',...
    'MaxGenerations',100,'MaxTime',60*60*5, 'PopulationSize', 50);
options.InitialPopulation = [
    
   13.0000   49.0000   34.0000    4.2804    0.9125
    8.0000   38.0000   39.0000    1.9306    0.4731
   13.0000   48.0000   35.0000    2.2072    0.8355
    8.0000   49.0000   34.0000    3.6982    0.9744
    6.0000   26.0000    9.0000    6.4558    0.6207
    7.0000   10.0000   14.0000    3.4775    0.6345
    8.0000   38.0000   39.0000    2.0711    0.4677
    8.0000   49.0000   35.0000    2.9931    0.7869
    8.0000   29.0000   38.0000    2.6834    0.6378
    7.0000   31.0000    6.0000    5.6905    0.7388
    7.0000   29.0000   14.0000    2.9309    0.5318
    6.0000   27.0000    9.0000    6.4589    0.6189
    1.0000    2.0000    1.0000    5.7644    0.7685
    6.0000    8.0000    3.0000    6.3160    0.9434
    8.0000   30.0000   38.0000    2.6887    0.5417
    8.0000   30.0000   38.0000    2.5937    0.6355
    8.0000   38.0000   39.0000    2.0236    0.4720
    2.0000    4.0000    1.0000    6.2784    0.9389


    ];
% Add initial guess into initial population matrix to speed up process
[x,fval,exitflag,output,population,scores] = ...
    gamultiobj(@COMPEL_2023_5Phase_BuckHybridBuck_D_Sweep_GAM, Nvars, A, b, Aeq, beq, LB.*sf, UB.*sf, [], INTCON, options);
% save('gaResults.mat', 'x', 'fval', 'exitflag', 'output','population','scores')
% load('gaResults.mat');
%}



%{
   26.0000   49.0000   30.0000    5.5938    1.0600   % This on has an issue Need to look into this
%}






