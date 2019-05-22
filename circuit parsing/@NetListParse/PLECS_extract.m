function [] = PLECS_extract(obj)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% PLECS Blockset input to AURA

% Ensure have the model file location within the MATLAB path
% First load the circuit

filename = obj.filename;

if ~(strcmp(filename(end-3,end),'.slx'))
  error('Did not recognize file extension for Simulink file');
end

filename_noext = filename(1:end-4);

load_system(filename);

% Assume path to circuit is filname/Circuit:
filename_Cir=strcat(filename_noext,'/Circuit');

names = plecs('get',filename_Cir,'StateSpaceOrder');

%% Get binary representation of number of states to change R and C for D and M
number_of_states = 2^length(names.Switches);
bin=de2bi(0:number_of_states-1);
state = bin;

% Iterate through each state and extract
for i = 1:1:number_of_states
    plecs('set',filename_Cir,'SwitchVector',state(i,:));
    t = plecs('get',filename_Cir,'Topology');
    Anum(:,:,i) = t.A;
    Bnum(:,:,i) = t.B;
    Cnum(:,:,i) = t.C;
    Dnum(:,:,i) = t.D;
    Inum(:,:,i) = t.I;
end

% Match extracted values to those in class
obj.Anum = Anum;
obj.Bnum = Bnum;
obj.Cnum = Cnum;
obj.Dnum = Dnum;

obj.Inum = Inum;

obj.InputNames = names.Inputs;
obj.StateNames = names.States;
obj.OutputNamesCD = names.Outputs;
obj.SwitchNames = names.Switches;

%%% Next step here is integrating these results into a simulation class and make adjustments to those functions %%%

% strcat(filename_noext,'.mat');
% save('Buck_Qual.mat');

%% How to perform a steady-state simulation in PLECS via MATLAB commands

% % For this example the simulink model is called 'plBuckSweep'
% % The PLECS circuit in the model is called 'Circuit'

% % Ensure have the model file location within the MATLAB path

% % Load the system using:
% load_system(plBuckSweep);

% % Then in order to conduct a steady-state simulation:
% % Note: The PLECS Steady-State Analysis simulation block must be within the model
% plsteadystate('plBuckSweep/Steady-State Analysis');

% % There a number of parameters that can be set in PLECS to affect how the
% % steady state simulation in conducted.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
% Parameter             Description
%
% TimeSpan              For a fixed system period, the period length;
%                       this is the least common multiple of the periods
%                       of independent sources in the system.
%                       For a variable system period, the maximum
%                       time span during which to look for a trigger
%                       event marking the end of a period. Set by the
%                       System period length/Max simulation time
%                       span field.
%
% TStart                Simulation start time. Set by the Simulation
%                       start time field.
%
%
% Tolerance             Relative error tolerance used in the convergence
%                       criterion. Set by the Termination tolerance
%                       field.
%
% MaxIter               Maximum number of iterations allowed. Set by
%                       the Max number of iterations field.
%
% Display               Specifies the level of detail of the diagnostic
%                       messages displayed in the command window
%                       (iteration, final, off). Set by the Display
%                       drop-down list.
%
% HideScopes            Hide all Simulink scope windows during an
%                       analysis in order to save time.
%
% HiddenStates          Specifies how to handle Simulink blocks with
%                       'hidden' states, i.e. states that are not stored in
%                       the state vector (error, warning, none). Set by
%                       the Hidden model states drop-down list.
%
% FinalStateName        Name of a MATLAB variable used to store the
%                       steady-state vector at the end of an analysis.
%                       Set by the Steady-state variable field.
%
% NCycles               Number of steady-state cycles that should be
%                       simulated at the end of an analysis. Set by the
%                       Show steady-state cycles field.
%
% JPert                 Relative perturbation of the state variables
%                       used to calculate the approximate Jacobian
%                       matrix.
%
% JacobianCalculation   Controls the way the Jacobian matrix is calculated
%                       (full, fast). The default is fast.
%
% NInitCycles           Number of cycle-by-cycle simulations that
%                       should be performed before the actual steady-state
%                       analysis. This parameter can be used to
%                       provide the algorithm with a better starting
%                       point. The default is 0.
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
