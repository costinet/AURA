%% TestExampleTopologies loads and simulates steady-state for a variety of example converters
%   The script requires properly-formatted netlist files for the example
%   topologies that include all initial parameters and inputs.  
%
%   if the allAssess variable is set, all of the topologies listed will be
%   simulated sequentially.  Otherwise, only those indices in the models 
%   varible will be simulated


clear all;

summaryStrings = {};

debug = 0;
allAssess = 1;      % run all models sequentially

models = 8;         % if allAssess == 0, this model will be run


%% Add test netlist folder to the path
sdir = mfilename('fullpath');
sdir = sdir(1:find(sdir=='\',1,'last')-1);
addpath(sdir);
addpath([sdir '\LTspice NetLists'])

%% Test Circuits
modelfile{1} = 'BuckNoDiodes.net'; 
modelfile{2} = 'BuckTest1.net'; 
modelfile{3} = 'Buck_Vout_D.net'; 
modelfile{4} = 'Buck_Vout_Dep.net'; 
modelfile{5} = 'Boost_Async.net'; 
modelfile{6} = 'Boost_Async_BD.net'; 
modelfile{7} = '3levelbuck.net'; 
modelfile{8} = '3levelbuckDep.net';
modelfile{9} = 'DAB.net'; 
modelfile{10} = 'DAB_secLeakage.net'; 
modelfile{11} = 'SeriesCapBuck.net'; 
modelfile{12} = 'WPT_VaryingU.net'; 
% modelfile{13} = 'Flyback.net'; 

modelfile(cellfun(@isempty,modelfile)) = [];


%% Most recent full run results:
% Model BuckNoDiodes.net converged after 6 iterations in 0.16597 seconds
% Model BuckTest1.net converged after 6 iterations in 0.073186 seconds
% Model Buck_Vout_D.net converged after 6 iterations in 0.034889 seconds
% Model Buck_Vout_Dep.net converged after 6 iterations in 0.013483 seconds
% Model Boost_Async.net converged after 30 iterations in 0.068542 seconds
% Model Boost_Async_BD.net converged after 9 iterations in 0.017218 seconds
% Model 3levelbuck.net converged after 5 iterations in 0.035796 seconds
% Model 3levelbuckDep.net converged after 5 iterations in 0.020215 seconds
% Model DAB.net converged after 7 iterations in 0.01748 seconds
% Model DAB_secLeakage.net converged after 8 iterations in 0.019213 seconds
% Model SeriesCapBuck.net converged after 6 iterations in 0.027859 seconds
% Model WPT_VaryingU.net converged after 8 iterations in 0.02606 seconds
% Model Flyback.net DID NOT converge within 102 in 10.4729 seconds

if allAssess
    % warning('off','all');
    h = waitbar(0,'Initializing');
    updateWaitBar(h, modelfile, 0, 0, 0, 0);
    models = 1:length(modelfile);
else 
    warning('on','all');
end

for selectedModel = models
        

    %% Load model file
    circuitPath = [modelfile{selectedModel}];

    
    %% Analyze circuit
    clear sim u us;
    sim = SMPSim(circuitPath);
    % [top, conv] = sim.initialize(circuitPath); 

    % If Debugging
        % sim.debug = 1;
        % sim.debug2 = 1;

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
        summaryStrings{selectedModel} = ['Model ' modelfile{selectedModel} ' converged after ' num2str(niter+1) ' iterations in ' num2str(solveTime) ' seconds'];
    elseif allAssess
        summaryStrings{selectedModel} = ['Model ' modelfile{selectedModel} ' DID NOT converge within ' num2str(niter+1) ' in ' num2str(solveTime) ' seconds'];
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
    
