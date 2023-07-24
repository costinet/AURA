%% TestExampleTopologies loads and simulates steady-state for a variety of example converters
%   The script requires propoerly-formatted simulink files for the example
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

models = 4;


if ~debug
    w = warning ('off','all');
else
    w = warning ('on','all');
end

sdir = mfilename('fullpath');
sdir = sdir(1:find(sdir=='\',1,'last')-1);
addpath(sdir);
addpath([sdir '\LTspice NetLists'])

%% Load test circuit
modelfile{1} = '3levelbuck.net'; 


%% Most recent full run results:



if allAssess
    warning('off','all');
    h = waitbar(0,'Initializing');
    updateWaitBar(h, modelfile, 0, 0, 0, 0);
    models = 1:length(modelfile);
else 
    warning('on','all');
end

for selectedModel = models
        
% %     % find_system(modelfile,'SearchDepth',1, 'IncludeCommented', 'on')
% %     open_system(modelfile{selectedModel},'loadonly');
    circuitPath = [modelfile{selectedModel}];
% %     set_param(circuitPath,'Commented','on');
% %     clear sim;
% %     simout = sim(modelfile{selectedModel},eps);
% %     
% %     for i = 1:length(simout.properties)
% %         assignin('base',simout.properties{i},eval(['simout.' simout.properties{i}]));
% %     end
% %     
% %     set_param(circuitPath,'Commented','off');
    
    
    %% Analyze circuit
    
    sim = SMPSim();
    conv = sim.converter;
    top = sim.topology;

%     sim.debug = 1;
%     sim.debug2 = 1;

swvec = [1 0 1 0; 0 1 0 1];
    
    top.loadCircuit(circuitPath,swvec,1);
    sim.u = us';
    conv.setSwitchingPattern(swvec, ts)

    if(debug)
        sim.steadyState;
        sim.plotAllStates(4)
    end

    
    niter = sim.findValidSteadyState;

    if(debug)
        sim.plotAllStates(1);
    end

    if allAssess
        updateWaitBar(h, modelfile, selectedModel, niter, 0, 0);
    end



    if(allAssess && niter <= 50)
        summaryStrings{selectedModel} = ['Model ' PLECsModel{selectedModel} ' converged after ' num2str(niter+1) ' iterations'];
    elseif allAssess
        summaryStrings{selectedModel} = ['Model ' PLECsModel{selectedModel} ' DID NOT converge within ' num2str(niter+1) ' iterations'];
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
    
