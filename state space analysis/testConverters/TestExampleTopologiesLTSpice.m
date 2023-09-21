%% TestExampleTopologies loads and simulates steady-state for a variety of example converters
%   The script requires properly-formatted netlist files for the example
%   topologies that include all initial parameters and inputs.  
%
%   if the allAssess variable is set, all of the topologies listed will be
%   simulated sequentially.  Otherwise, only those indices in the models 
%   varible will be simulated


clear all;

summaryStrings = {};

debug = 1;
allAssess = 1;      % run all models sequentially

models = 6;         % if allAssess == 0, this model will be run


%% Add test netlist folder to the path
sdir = mfilename('fullpath');
sdir = sdir(1:find(sdir=='\',1,'last')-1);
addpath(sdir);
addpath([sdir '\LTspice NetLists'])

%% Test Circuits
modelfile{1} = 'BuckNoDiodes.net'; 
modelfile{2} = 'BuckTest1.net'; 
modelfile{3} = 'Buck_Vout_D.net'; 
modelfile{4} = '3levelbuck.net'; 
modelfile{5} = '3levelbuckDep.net'; 
modelfile{6} = 'DAB.net'; 


%% Most recent full run results:
% Model BuckNoDiodes.net converged after 15 iterations in 0.34148 seconds
% Model BuckTest1.net converged after 14 iterations in 0.23347 seconds
% Model Buck_Vout_D.net converged after 11 iterations in 0.10838 seconds
% Model 3levelbuck.net converged after 8 iterations in 0.10606 seconds
% Model 3levelbuckDep.net converged after 8 iterations in 0.087727 seconds
% Model DAB.net converged after 24 iterations in 0.56207 seconds


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
    circuitPath = [modelfile{selectedModel}];

    
    %% Analyze circuit
    sim = SMPSim();
    [top, conv] = sim.initialize(circuitPath);   

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
    
