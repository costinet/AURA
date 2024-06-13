%% TestExampleTopologies loads and simulates steady-state for a variety of example converters
%   The script requires properly-formatted simulink files for the example
%   topologies that include all initial parameters and inputs.  These are
%   loaded into the base workspace so that PLECS will see them when the
%   state space matrices are parsed.
%
%   if the allAssess variable is set, all of the topologies listed will be
%   simulated sequentially.  Otherwise, only those indices in the models 
%   varible will be simulated


clear all;

summaryStrings = {};

debug = 0;
allAssess = 1;      % run all models sequentially

models = 13;         % if allAssess == 0, this model will be run

%% Add test circuits folder to the path
sdir = mfilename('fullpath');
sdir = sdir(1:find(sdir=='\',1,'last')-1);
addpath(sdir);
addpath([sdir '\ExampleTopologies'])

%% Test Circuits
modelfile{1} = 'AsyncBoost'; PLECsModel{1} = 'Boost_Async';
modelfile{2} = 'MRBuckDiodes'; PLECsModel{2} = 'MRBuck';
modelfile{3} = 'DSC4to1'; PLECsModel{3} = 'HDSC';
modelfile{4} = 'DAB'; PLECsModel{4} = 'DAB_oneCap';
modelfile{5} = 'DABfull'; PLECsModel{5} = 'DAB_8Cap';
modelfile{6} = 'DSC4to1Diodes'; PLECsModel{6} = 'HDSC_withDiodes';
modelfile{7} = 'SCBuckHybridBuck'; PLECsModel{7} = 'BuckBuck';
modelfile{8} = 'SCFibb2S_diodes'; PLECsModel{8} = 'Q-FibonacciESRCossDiodes';
modelfile{9} = 'SCFibb1S_diodes'; PLECsModel{9} = 'Q-FibonacciESRCossDiodes';
modelfile{10} = 'TwoPhaseBuck'; PLECsModel{10} = 'TwoAsymmBuck';
modelfile{11} = 'SCBuckHybridBuck_COMPEL23'; PLECsModel{11} = 'BuckBuck';
modelfile{12} = 'DAB_R'; PLECsModel{12} = 'DAB_Rload';
modelfile{13} = 'WPT_varyingU'; PLECsModel{13} = 'WPT_SSP';


% modelfile{length(modelfile)+1} = 'DABFullNoDiodes'; PLECsModel{length(PLECsModel)+1} = 'DAB_8Cap';
% modelfile{length(modelfile)+1} = 'SCFibb2S'; PLECsModel{length(PLECsModel)+1} = 'Q-FibonacciESRCoss';
% modelfile{length(modelfile)+1} = 'Flyback'; PLECsModel{length(PLECsModel)+1} = 'PC_Flyback';

%% Most recent full run results:
% Model Boost_Async converged after 43 iterations in 0.39162 seconds
% Model MRBuck converged after 11 iterations in 0.088323 seconds
% Model HDSC converged after 6 iterations in 0.13764 seconds
% Model DAB_oneCap converged after 37 iterations in 0.74975 seconds
% Model DAB_8Cap DID NOT converge within 102 in 1.7696 seconds
% Model HDSC_withDiodes converged after 2 iterations in 0.0043238 seconds
% Model BuckBuck converged after 9 iterations in 0.79623 seconds
% Model Q-FibonacciESRCossDiodes converged after 8 iterations in 0.48962 seconds
% Model Q-FibonacciESRCossDiodes converged after 2 iterations in 0.0027774 seconds
% Model TwoAsymmBuck converged after 24 iterations in 0.222 seconds
% Model BuckBuck converged after 15 iterations in 5.0886 seconds
% Model DAB_Rload converged after 41 iterations in 0.931 seconds
% Model WPT_SSP converged after 11 iterations in 0.14177 seconds
%
% Model DAB_8Cap DID NOT converge within 102 iterations
% Model Q-FibonacciESRCoss converged after 8 iterations



if allAssess
    warning('off','all');
    h = waitbar(0,'Initializing');
    updateWaitBar(h, modelfile, 0, 0, 0, 0);
    models = 1:length(modelfile);
else 
    warning('on','all');
end

for selectedModel = models
        
    % If Debugging
    %     sim.debug = 1;
    %     sim.debug2 = 1;

    %% Load model file
    circuitPath = [modelfile{selectedModel} '/' PLECsModel{selectedModel}];
    loadSimulinkOutputsToBase(modelfile{selectedModel},PLECsModel{selectedModel});
    
    
    %% Analyze circuit
    sim = SMPSim();
    sim.initialize(circuitPath, swvec, us, ts);

    if(debug)
        sim.steadyState;
        sim.plotAllStates(4)
    end

    tic;
    niter = sim.findValidSteadyState;   
    solveTime = toc;

    if(debug)
        sim.plotAllStates(1);
    end

    %% Update progress bar and store result

    if allAssess
        updateWaitBar(h, modelfile, selectedModel, niter, 0, 0);
    end

    if(allAssess && niter <= sim.maxItns)
        summaryStrings{selectedModel} = ['Model ' PLECsModel{selectedModel} ' converged after ' num2str(niter+1) ' iterations in ' num2str(solveTime) ' seconds'];
    elseif allAssess
        summaryStrings{selectedModel} = ['Model ' PLECsModel{selectedModel} ' DID NOT converge within ' num2str(niter+1) ' in ' num2str(solveTime) ' seconds'];
    end
end

for i = 1:length(summaryStrings)
    disp(summaryStrings{i})
end

if allAssess
    close(h);
end

function updateWaitBar(h, modelfile, selectedModel, niter, altered, finalRun)
    message = modelfile;
    for i = 1:length(modelfile)
        if i < selectedModel
            message{i} = [modelfile{i} '... completed'];
        elseif i > selectedModel
            message{i} = [modelfile{i} '... waiting'];
        else
            message{i} = [modelfile{i} '... running.  Iteration ' num2str(niter) '/250, altered = ' num2str(altered) ' finalRun = ' num2str(finalRun)];
        end
    end
    waitbar(selectedModel/length(modelfile),h,message)
    h.Position(4) = length(modelfile)*18;
end
    
