clearvars;

%% Add test circuits folder to the path
sdir = mfilename('fullpath');
sdir = sdir(1:find(sdir=='\',1,'last')-1);
addpath(sdir);
addpath([sdir '\ExampleTopologies'])

%% Example Steady-State Solve from PLECS circuit
modelfile = 'AsyncBoost'; PLECsModel = 'Boost_Async';

loadSimulinkOutputsToBase(modelfile,PLECsModel);
circuitPath = [modelfile '/' PLECsModel];

sim = SMPSim();
sim.initialize(circuitPath);                    %swvec, us, and ts are defined 
                                                %in the simulink file and
                                                %loaded via loadSimulinkOutputsToBase
sim.steadyState();
sim.plotAllStates(1)
sim.plotAllOutputs(2)

%% Example Steady-State Solve from LTspice netlist
modelfile = 'BuckTest1.net';

sim = SMPSim();
sim.initialize(modelfile);                      %swvec, us, and ts are defined 
                                                %in the netlist file and
                                                %loaded during parsing
sim.steadyState();
sim.plotAllStates(3)
sim.plotAllOutputs(4)

%% Example Steady-State Solve from PLECS circuit with dependent States
modelfile = 'AsyncBoost'; PLECsModel = 'Boost_Async';

loadSimulinkOutputsToBase(modelfile,PLECsModel);
circuitPath = [modelfile '/' PLECsModel];

sim = SMPSim();
sim.initialize(circuitPath, swvec, us, ts);     %here swvec, us, and ts are explicitly specified

sim.suppressIterationOutput = 1;    %Suppress per-iteration output to command line

sim.plotAllStates(1);               %Plot all states BEFORE solving valid diode conduction
hold(gcf().Children,'on');          %hold on, to compare with final solution

sim.findValidSteadyState();         %solves correct dependent state behaviour
sim.plotAllStates(1)                %Plot all states AFTER solving valid diode conduction



%% Example Steady-State Solve from LTspice netlist with dependent States
modelfile = 'BuckTest1.net';

sim = SMPSim();
sim.initialize(modelfile);                      %swvec, us, and ts are defined 
                                                %in the netlist file and
                                                %loaded during parsing
sim.suppressIterationOutput = 1;    %Suppress per-iteration output to command line

sim.plotAllStates(3);               %Plot all states BEFORE solving valid diode conduction
hold(gcf().Children,'on');          %hold on, to compare with final solution

sim.findValidSteadyState();         %solves correct dependent state behaviour
sim.plotAllStates(3)                %Plot all states AFTER solving valid diode conduction

%% Example Underdefined PLECS File
modelfile = 'DAB_Incomplete'; PLECsModel = 'DAB_oneCap';

open_system(modelfile,'loadonly');
circuitPath = [modelfile '/' PLECsModel];
    
sim = SMPSim();
conv = sim.converter;
top = sim.topology;
    
try
    top.loadCircuit(circuitPath); %This fails with explanation of missing variables
catch e
    disp(getReport( e, 'extended', 'hyperlinks', 'on' ) )
end

