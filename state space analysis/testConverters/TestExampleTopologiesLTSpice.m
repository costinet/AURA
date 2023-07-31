%% TestExampleTopologies loads and simulates steady-state for a variety of example converters
%   The script requires properly-formatted netlist files for the example
%   topologies that include all initial parameters and inputs.  
%
%   if the allAssess variable is set, all of the topologies listed will be
%   simulated sequentially.  Otherwise, only those indices in the models 
%   varible will be simulated


clear all;

summaryStrings = {};

allAssess = 1;      % run all models sequentially

models = 1;         % if allAssess == 0, this model will be run


%% Add test netlist folder to the path
sdir = mfilename('fullpath');
sdir = sdir(1:find(sdir=='\',1,'last')-1);
addpath(sdir);
addpath([sdir '\LTspice NetLists'])

%% Test Circuits
modelfile{1} = 'BuckTest1.net'; 
modelfile{2} = 'Buck_Vout_D.net'; 


%% Most recent full run results:
% Model BuckTest1.net converged after 18 iterations in 0.53823 seconds
% Model Buck_Vout_D.net converged after 11 iterations in 0.15198 seconds


if allAssess
    warning('off','all');
    h = waitbar(0,'Initializing');
    updateWaitBar(h, modelfile, 0, 0, 0, 0);
    models = 1:length(modelfile);
else 
    warning('on','all');
end

for selectedModel = models
        

    circuitPath = [modelfile{selectedModel}];

    
    %% Analyze circuit
    
    sim = SMPSim();
    conv = sim.converter;
    top = sim.topology;

%     sim.debug = 1;
%     sim.debug2 = 1;


    
    top.loadCircuit(circuitPath,[],1);
    sim.u = us;
    conv.setSwitchingPattern(swvec, ts)

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
    
