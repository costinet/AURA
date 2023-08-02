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
UB = [ 80 80 80 50    5e6];
sf = [  1  1  1  0.1  1e-6];  %scale factor to ~ equalize magnitudes of variables

Aeq = [];
beq = [];
A = [];
b = [];

INTCON = [1 2 3];


%% Genetic Algorithm

options = optimoptions('gamultiobj','PlotFcn',{@gaplotpareto,@gaplotparetodistance,@gaplotrankhist,@gaplotspread},...
    'TolCon',1e-9,'Display','iter',...
    'MaxGenerations',75,'MaxTime',60*60*6, 'PopulationSize', 75);
options.InitialPopulation = [
    6.0000   21.0000   25.0000    1.8995    0.6721
    6.0000   10.0000    1.0000    1.6763    1.1747
    6.0000   21.0000   29.0000    1.7739    0.5753
    6.0000    1.0000   29.0000    1.5295    0.8948
    8.0000   21.0000   37.0000    1.2560    0.6483
    6.0000   21.0000    5.0000    1.3150    0.6422
    7.0000   22.0000   58.0000    1.0965    0.3840
    6.0000    1.0000    1.0000    2.0000    0.7457
    6.0000    8.0000    1.0000    1.9268    0.7511
    8.0000   21.0000   34.0000    1.1875    0.4450
    8.0000   21.0000   34.0000    1.2699    0.4904
    6.0000    1.0000    1.0000    1.7427    1.0407
    6.0000   21.0000   34.0000    1.2782    0.5380
    7.0000   21.0000   33.0000    1.7291    0.5585
    8.0000   21.0000   37.0000    1.1332    0.5061
    8.0000   22.0000   37.0000    1.3723    0.5004
    8.0000   22.0000   37.0000    1.3121    0.4976
    6.0000    1.0000    1.0000    1.5587    0.8222
    7.0000   22.0000   58.0000    1.0776    0.3676
    6.0000    6.0000    1.0000    1.9934    0.5819
    6.0000    7.0000    1.0000    1.9286    0.7917
    7.0000   21.0000   33.0000    1.5751    0.5587
    8.0000   22.0000   37.0000    1.4468    0.5135
    6.0000    7.0000    1.0000    1.9877    0.8513
    7.0000   21.0000   34.0000    1.3034    0.5464
    7.0000   21.0000   37.0000    1.1331    0.4860
    8.0000   22.0000   37.0000    1.2812    0.5011


    ];
% Add initial guess into initial population matrix to speed up process
[x,fval,exitflag,output,population,scores] = ...
    gamultiobj(@COMPEL_2023_HybridDickson_D_Sweep_GAM, Nvars, A, b, Aeq, beq, LB.*sf, UB.*sf, [], INTCON, options);
% save('gaResults.mat', 'x', 'fval', 'exitflag', 'output','population','scores')
% load('gaResults.mat');
%}



%{
   26.0000   49.0000   30.0000    5.5938    1.0600   % This on has an issue Need to look into this
%}






