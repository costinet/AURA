

clear

%%{
%% Set up optimization for 9 Variables
Nvars = 10;
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


LB = [1 1 1 1 1 1 1 0.5e6 0.0007 0.0007];
UB = [11 11 11 11 11 11 4 2e6 0.008 0.008 ];
sf = [1 1 1 1 1 1 1 1e-6 100 100];  %scale factor to ~ equalize magnitudes of variables




Aeq = [];
beq = [];
A = [];
b = [];
INTCON = [1 2 3 4 5 6 7];

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
    'MaxGenerations',500,'MaxTime',60*60*2.5, 'MaxStallGenerations',50, 'PopulationSize', 50);
options.InitialPopulation = [
    
4.0000    4.0000    4.0000    4.0000    4.0000    4.0000  1   1.0e6       0.005       0.005
9.0000    4.0000    6.0000    6.0000    4.0000    4.0000  1  1.3078e6    0.004211    0.003431
4.0000    4.0000    4.0000    4.0000    2.0000    2.0000  1  1.4929e6    0.002722    0.006690
4.0000    4.0000    4.0000    4.0000    2.0000    2.0000  1  1.8182e6    0.002322    0.006922
4.0000    4.0000    4.0000    4.0000    2.0000    2.0000  1  1.8166e6    0.002305    0.007238
4.0000    4.0000    4.0000    4.0000    2.0000    2.0000  1  1.8043e6    0.002294    0.007498
4.0000    4.0000    4.0000    4.0000    2.0000    2.0000  1  1.8108e6    0.002297    0.007531
4.0000    4.0000    4.0000    4.0000    2.0000    2.0000  1  1.8167e6    0.002322    0.007576
4.0000    2.0000    2.0000    2.0000    2.0000    2.0000  1  1.8223e6    0.003260    0.007341
2.0000    2.0000    2.0000    2.0000    2.0000    2.0000  1  1.9162e6    0.003268    0.000887
2.0000    2.0000    2.0000    2.0000    2.0000    2.0000  1  1.6301e6    0.006522    0.007939
4.0000    4.0000    4.0000    4.0000    2.0000    2.0000  1  1.8175e6    0.002348    0.007599
2.0000    2.0000    2.0000    2.0000    2.0000    2.0000  1  1.6854e6    0.007417    0.008000
].*sf;
% Add initial guess into initial population matrix to speed up process
[x,fval,exitflag,output,population,scores] = ...
    ga(@AURA_Eff_Sweep_Disc_Buck_3L_Boost, Nvars, A, b, Aeq, beq, LB.*sf, UB.*sf, [], INTCON, options);
% save('gaResults.mat', 'x', 'fval', 'exitflag', 'output','population','scores')
% load('gaResults.mat');
%}



%{

3 level buck boost

4.0000    4.0000    4.0000    4.0000    4.0000    4.0000    1.0e6       0.005       0.005
9.0000    4.0000    6.0000    6.0000    4.0000    4.0000    1.3078e6    0.004211    0.003431
4.0000    4.0000    4.0000    4.0000    2.0000    2.0000    1.4929e6    0.002722    0.006690
4.0000    4.0000    4.0000    4.0000    2.0000    2.0000    1.8182e6    0.002322    0.006922
4.0000    4.0000    4.0000    4.0000    2.0000    2.0000    1.8166e6    0.002305    0.007238
4.0000    4.0000    4.0000    4.0000    2.0000    2.0000    1.8043e6    0.002294    0.007498
4.0000    4.0000    4.0000    4.0000    2.0000    2.0000    1.8108e6    0.002297    0.007531
4.0000    4.0000    4.0000    4.0000    2.0000    2.0000    1.8167e6    0.002322    0.007576
4.0000    2.0000    2.0000    2.0000    2.0000    2.0000    1.8223e6    0.003260    0.007341
2.0000    2.0000    2.0000    2.0000    2.0000    2.0000    1.9162e6    0.003268    0.000887
2.0000    2.0000    2.0000    2.0000    2.0000    2.0000    1.6301e6    0.006522    0.007939
4.0000    4.0000    4.0000    4.0000    2.0000    2.0000    1.8175e6    0.002348    0.007599
2.0000    2.0000    2.0000    2.0000    2.0000    2.0000    1.6854e6    0.007417    0.008000


%}

L_values=[100e-9 125e-9 150e-9 175e-9 200e-9 225e-9 250e-9 275e-9 300e-9 500e-9 750e-9 1e-6 2.5e-6 5e-6 7.5e-6 10e-6];


sf = [1 1 1 1 1 1 1e6 1e-6 100 100];
X = [4.0000    4.0000    4.0000    4.0000    2.0000    2.0000  100e-9   1.8108e6    0.002297    0.007531];

for i = 1:length(L_values)
    X(7) = L_values(i);
    Fval_L_sweep(i) = AURA_Eff_Sweep_Disc_Buck_3L_Boost(X.*sf);
end

