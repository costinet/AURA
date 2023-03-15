clear all;

summaryStrings = {};

debug = 0;
% % debug2 = 0;
% % finalRunMethod = 0; %1 to use timestepping, 0 to use eigs2tis
% % timeSteppingInit = 0;
allAssess = 1;      % run all models sequentially

% tmult = 0.5;

if ~debug
    w = warning ('off','all');
else
    w = warning ('on','all');
end

sdir = mfilename('fullpath');
sdir = sdir(1:find(sdir=='\',1,'last')-1);
addpath(sdir);

%% Load test circuit
modelfile{1} = 'AsyncBoost'; PLECsModel{1} = 'Boost_Async';
modelfile{2} = 'MRBuck'; PLECsModel{2} = 'MRBuck';
modelfile{3} = 'DSC4to1'; PLECsModel{3} = 'HDSC';
modelfile{4} = 'DAB'; PLECsModel{4} = 'DAB_oneCap';
modelfile{5} = 'DABfull'; PLECsModel{5} = 'DAB_8Cap';
modelfile{6} = 'DSC4to1Diodes'; PLECsModel{6} = 'HDSC_withDiodes';
modelfile{7} = 'Flyback'; PLECsModel{7} = 'PC_Flyback';
modelfile{8} = 'SCBuckHybridBuck'; PLECsModel{8} = 'BuckBuck';
modelfile{9} = 'DABFullNoDiodes'; PLECsModel{9} = 'DAB_8Cap';
modelfile{10} = 'SCFibb2S'; PLECsModel{10} = 'Q-FibonacciESRCoss';
modelfile{11} = 'SCFibb2S_diodes'; PLECsModel{11} = 'Q-FibonacciESRCossDiodes';
modelfile{12} = 'SCFibb1S_diodes'; PLECsModel{12} = 'Q-FibonacciESRCossDiodes';

models = 8;

if allAssess
    warning('off','all');
    h = waitbar(0,'Initializing');
    updateWaitBar(h, modelfile, 0, 0, 0, 0);
    models = 1:length(modelfile);
else 
    warning('on','all');
end

for selectedModel = models
        
    % find_system(modelfile,'SearchDepth',1, 'IncludeCommented', 'on')
    open_system(modelfile{selectedModel},'loadonly');
    circuitPath = [modelfile{selectedModel} '/' PLECsModel{selectedModel}];
    set_param(circuitPath,'Commented','on');
    clear sim;
    simout = sim(modelfile{selectedModel},eps);
    
    for i = 1:length(simout.properties)
        assignin('base',simout.properties{i},eval(['simout.' simout.properties{i}]));
    end
    
    set_param(circuitPath,'Commented','off');
    
    
    %% Analyze circuit
    
    sim = SMPSim();
    conv = sim.converter;
    top = sim.topology;
    
    top.loadCircuit(circuitPath,swvec,1);
    sim.u = us';
    conv.setSwitchingPattern(1:size(swvec,1), ts)
    
    sim.steadyState;
    if(debug)
        sim.plotAllStates(1);
    end
    
    niter = sim.findValidSteadyState;

    updateWaitBar(h, modelfile, selectedModel, niter, 0, 0);

% ssOrder = plecs('get', circuitPath, 'StateSpaceOrder');
% outputs = ssOrder.Outputs;

% % outputs = top.outputLabels;

% % %% finalRun 
% % % once everything seems to be error-free based on discrete time points,
% % % goes through once more with eigenvalue-based spacing to make sure no 
% % % inter-sample violations are occuring.  
% % finalRun = 0;
% % 
% % %% Symmetry check
% % % May be useful but not doing anything with it yet.  Can identify that DAB,
% % % etc. exhibit half-cycle symmetry
% % TF = conv.checkForSymmetry;
% % 
% % %%  TimeStepping Attempt
% % if timeSteppingInit > 0
% %     tic
% %     [Xf,ts,swinds] = timeSteppingPeriod(sim);
% %     
% %     numPeriods = timeSteppingInit;
% %     Xsims = zeros(size(Xf,1),numPeriods);
% %     for i = 1:numPeriods
% %         conv.setSwitchingPattern(swinds, ts);
% %     %     Xss = sim.steadyState;
% %         Xsims(:,i) = Xf;
% %         [Xf,ts,swinds] = sim.timeSteppingPeriod(Xf, ts, swinds );
% %         if debug == 1
% %             disp([Xsims(:,i) Xf]);
% %         end
% %     end
% % 
% %     [newts,newswinds] = sim.format1DtimingVector(ts,swinds);
% %     
% % 
% % %     switchedDevs = diff(top.swseq(swinds,:),1);
% % %     diodes = ~cellfun(@isempty,regexp(top.switchLabels, 'D.*'));
% % %     activeSwitched = switchedDevs(:,~diodes);
% % %     newInt = sum(abs(activeSwitched),2)~=0;
% % % 
% % %     newts = zeros(length(newInt)+1, sum(newInt)+1);
% % %     newswinds = zeros(length(newInt)+1, sum(newInt)+1);
% % % 
% % %     newts(1) = ts(1);
% % %     newswinds(1) = swinds(1);
% % %     col = 1;
% % %     row = 2;
% % %     for i = 1:length(newInt)
% % %         if newInt(i) == 1
% % %             col = col + 1;
% % %             row = 1;
% % %         end
% % %         newts(row,col) = ts(i+1);
% % %         newswinds(row,col) = swinds(i+1);
% % %         row = row+1;
% % %     end
% % % 
% % %     newts(all(newswinds == 0,2),:)=[];
% % %     newswinds(all(newswinds == 0,2),:)=[];
% % 
% %     
% % %     warning('This is not correct, generally');
% % %     ts = [[ts(1); zeros(length(ts)-2,1)], [ts(2:end)]];
% % %     swinds = [[swinds(1); zeros(length(swinds)-2,1)], [swinds(2:end)']];
% % 
% % %     load('flyBackTimeStepping')
% %     conv.setSwitchingPattern(newswinds, newts)
% %     clear newts;
% %     toc
% %     X0 = sim.Xs(:,1);
% %     % save('flyBackTimeStepping', 'Xf', 'swinds', 'swvec', 'ts')
% % end
% % 
% % tic
% % while(1)
% %     Xss = sim.steadyState;
% %     
% %     %% Update constraints per the current switching vector
% %     Cbnd = top.Cbnd; Dbnd = top.Dbnd; 
% %     hyst = top.bndHyst; switchRef = top.switchRef;
% %     Cbnd = Cbnd(:,:,conv.swind);
% %     Dbnd = Dbnd(:,:,conv.swind);
% % 
% %     if(0)
% %         %% ~~Continuous time violation margin
% %         [ xs, t] = sim.SS_WF_Reconstruct;
% %         for i = 1:length(xs)-1
% %             swstate = find(t(i) <= cumsum(ts),1,'first');
% %             violationMargin(:,i) = Cbnd(:,:,swstate)*xs(:,i) + Dbnd(:,:,swstate)*us' ...
% %                 - hyst(:,1) + hyst(:,2);
% %         end
% % 
% %         figure; plot(t(1:end-1),violationMargin)
% %         ylim([-100 100])
% %         legend(sim.switchNames) %% <- problem here, now.
% %     end
% % 
% %     %% Discrete timepoint violation margin
% %     [violateMarginStart,violateMarginEnd,targetValStart,targetValEnd] = sim.checkDiscreteErr;
% %     errBefore = min(violateMarginStart,0);
% %     errAfter = min(violateMarginEnd,0);
% % 
% % 
% %     if debug2 == 1
% %         sim.describeDiscreteErrors;
% %     end
% % 
% %     %% Check if we have the correct state to change times
% %     % if yes, proceed to jacobian.
% %     % if no, add states, then redo steady-state & error detection.
% %     % UPDATE:  just do it anyway, let addUncontrolledSwitching report on whether
% %     % it actually did anything.
% % 
% %     %% Insert additional states as necessary
% %     % We only need to insert a new state if 
% %     %   -- The interface has an error both before and after the switching
% %     %   -- OR it has one of the above, and the before and after switching
% %     %   positions aren't part of the modifiable uncontrolled times.
% %     
% %     [~,ints,~] = getIntervalts(conv);
% %     ints = ints';
% %     
% %     [tLocs,insertAt,adjType] = sim.findRequiredUncontrolledSwitching(violateMarginStart,violateMarginEnd);
% % 
% %     altered = 0;
% %     allChanges = [];
% % %     warning("Need something so it doesn't do corrections both before and after");
% % 
% %     for i = flip(find(insertAt))
% %         % addUncontrolledSwitching(obj, interval, beforeAfter, initialTime, switches, newStates)
% % %         dt = min( min(conv.controlledts)/20 , conv.ts(i)/20);
% %         [~, dts] = conv.getDeltaT();
% %         dt = min([min(dts),min(conv.controlledts)/20 , conv.ts(i)/20]);
% %         for j = 1:2
% %             if(any(adjType(:,i,j)))
% %                 [alt, newSwInd] = ...
% %                     conv.addUncontrolledSwitching(i,(-1)^(j+1), ...
% %                     dt,switchRef(tLocs(:,i),1),~switchRef(tLocs(:,i),2));
% %                 altered = altered | alt;
% % 
% %                 if debug2 == 1
% % %                     interval, beforeAfter, newSwInd, switches, newStates
% %                     locRefVec = [find(tLocs(:,i))]';
% %                     for locRef = locRefVec
% %                         allChanges = [allChanges; i, (-1)^(j+1), newSwInd, switchRef(locRef,1), ~switchRef(locRef,2)];
% %                     end
% %                     if i == find(insertAt,1,'first')
% %                         sim.describeInsertedIntervals(allChanges)
% %                     end
% %                 end
% %             end
% %         end
% %     end
% % 
% % 
% %     Xss = sim.steadyState;
% %     if(debug)
% %         sim.plotAllStates(1);
% %     end
% % 
% %     %% If anything was altered
% %     % in the above, we insert new switching sequences wherever needed to 
% %     % make sure we can (possibly) adjust timings to get zero error.  If we
% %     % changed the switching sequence, just loop back around until we have a
% %     % valid sequence.  Once the sequence isn't modified, proceed to adjust
% %     % timings.
% %     if ~altered
% %         if ~any(any(errBefore) | any(errAfter)) && finalRun == 0
% %             % If no errors in discrete time state vectors, run one more
% %             % time with higher time resolution to make sure no errors occur
% %             % in between discrete samples.
% %             finalRun = 1;
% %             if finalRunMethod
% %                 Xss = sim.steadyState;
% %                 [Xf,ts,swinds] = sim.timeSteppingPeriod();
% %                 conv.setSwitchingPattern(swinds, ts);
% %                 Xss = sim.steadyState;
% %             else
% %                 eigs2tis(conv);
% %                 Xss = sim.steadyState;
% %             end
% %             continue;
% %         elseif ~any(any(errBefore) | any(errAfter)) && finalRun == 1
% %             % If no errors and we've done the above, we're finished
% %             break;
% %         else
% %             % Otherwise keep looping
% %             finalRun = 0;
% %         end
% %         
% %         %% Correct based on discrete jacobian
% %         [JoutStart,JoutEnd] = sim.discreteJacobianConstraint;
% %         
% %         [i1,j1] = find(errBefore);
% %         [i2,j2] = find(errAfter);
% %         intV = [j1; j2]; 
% %         stateV = [i1; i2];
% % 
% %         A = zeros(length(j1)+length(j2),length(conv.ts));
% %         b = zeros(length(j1)+length(j2),1);
% %         tindex = 1:length(conv.ts);
% % 
% %         for i = 1:length(j1)
% %             deltas = squeeze(JoutStart(i1(i),j1(i),:));
% %             A(i,:) = deltas';
% %             b(i,1) = targetValStart(i1(i),j1(i));
% %         end
% %         if(isempty(i)), i=0; end
% %         for j = 1:length(j2)
% %             deltas = squeeze(JoutEnd(i2(j),j2(j),:));
% %             A(i+j,:) = deltas';
% %             b(i+j,1) = targetValEnd(i2(j),j2(j));
% %         end
% %         if(isempty(j)), j=0; end
% % 
% %         unChangeable = isnan(sum(A,1));
% % 
% %         %% Check for slope change (Nonfunctional)
% % %         interinterval = (errAfter~= 0 & errBefore == 0);
% % %         sign(JoutStart) ~= sign(JoutEnd);
% % 
% %        
% % 
% %         %% Attempt: add zero net perturbation to time as a part of the
% %         %% equations -- MAY NOT  because some times are dropped.
% %         scaleF = norm(A(:,~unChangeable))/numel(A);
% %         [~, timeInts, ~] = conv.getIntervalts;
% %         A = [A; zeros(max(timeInts), size(A,2))];
% %         b = [b; zeros(max(timeInts), 1)];
% %         for i = 1:max(timeInts)
% %             if ~any(timeInts'==i & unChangeable)
% %                 A(end-max(timeInts)+i,timeInts==i) = scaleF;
% %             end
% %         end
% %         
% %         A(:,unChangeable) = [];
% % 
% %         tsolve = zeros(size(conv.ts));
% %         tsolve(~unChangeable) = -(A\b);
% %         
% %         if any(isnan(tsolve))
% %             tsolve(~unChangeable) = -pinv(A)*b;
% %             warning('Sometimes this goes awry')
% % 
% % %             tsolve = tsolve*min(conv.ts)/max(abs(tsolve))*.01;
% %             
% %         end
% % 
% %         oldts = conv.ts;
% %         tps = conv.validateTimePerturbations2(tsolve);
% %         for i=length(tsolve):-1:1
% %             if tsolve(i) ~= 0
% %                 conv.adjustUncontrolledTiming(i, tps(i));
% %             end
% %         end
% % %         
% % %         [~,dtLims] = getDeltaT(sim.converter);
% % %         
% % %         tr = tsolve./dtLims;
% % %         tr = max(tr(tr>1));
% % %         
% % %         if ~isempty(tr)
% % %             tsolve = tsolve/tr;
% % %         end
% % %         
% %         
% % 
% %         
% % %         recurs = 0;
% % %         for i=length(tsolve):-1:1
% % %            if tsolve(i) ~= 0
% % % %                conv.adjustTiming(i, deltaTs(i));
% % %                 
% % %                 %% If adjustUncontrolledTiming has an issue and changes the interval
% % %                 % should use validateTimePerturbations to find this on the
% % %                 % front end
% % %                 if recurs ~=1
% % %                     [~, recurs] = conv.adjustUncontrolledTiming(i, tsolve(i));
% % %                 else
% % %                     recurs = 0;
% % %                 end
% % %            end
% % %         end
% % % % %         conv.adjustUncontrolledTimingVector(1:length(tsolve), tsolve)
% %         if exist('newts','var')
% %             if numel(newts) == numel(conv.ts)
% %                 if max(abs(newts - conv.ts)) < 10*conv.timingThreshold
% %                     % If the code gets in here, it looks like we're
% %                     % oscillating, so try to break out
% %                     [Xf,ts,swinds] = timeSteppingPeriod(sim);
% %                     [newts,newswinds] = sim.format1DtimingVector(ts,swinds);
% %                     conv.setSwitchingPattern(newswinds, newts)
% %                     warning('hack')
% %                 end
% %             end
% %         end
% % 
% %         newts = conv.ts;
% % 
% % %         [[newts-oldts]', tsolve']
% % 
% % 
% % 
% % %         warning('Can get stuck trying to make intervals longer when it cannot');
% % 
% %         
% %         if(debug)
% % %             disp(tsolve)
% % %             disp(sum(errBefore + errAfter, 'all'))
% %             Xss = sim.steadyState;
% %             sim.plotAllStates(10);
% %         end
% %         
% %         niter = niter+1;
% %         
% %         if(~allAssess)
% % %             weightTotalErr = getWeightedTotalError(sim, errBefore,errAfter);
% % %             disp([niter weightTotalErr]);
% %             disp([niter sum(errBefore + errAfter, 'all')]);
% %         else
% %             updateWaitBar(h, modelfile, selectedModel, niter, altered, finalRun);
% %         end
% %         
% %         if(allAssess && niter > 250)
% %             break;
% %         end
% %         
% %         if(~any(tsolve))
% %             error('timing not modified');
% %             break;
% %         end
% %         if(0)
% %             [ avgXs, avgYs ] = sim.ssAvgs(Xss);
% %             Pout = avgYs(21)*Vb1 + avgYs(22)*Vb2
% %             Pin = Vg*avgXs(3)
% %         end
% %         
% % 
% %     end
% % end
% % toc

    if(allAssess && niter <= 250)
        summaryStrings{selectedModel} = ['Model ' PLECsModel{selectedModel} ' converged after ' num2str(niter+1) ' iterations'];
    elseif allAssess
        summaryStrings{selectedModel} = ['Model ' PLECsModel{selectedModel} ' DID NOT converge within ' num2str(niter+1) ' iterations'];
    end
end

for i = 1:length(summaryStrings)
    disp(summaryStrings{i})
end
Xss = sim.steadyState;
sim.plotAllStates(10);
% sim.plotAllStates(11, [1, 15, 20, 11, 12, 7, 16, 21, 2]  );
% sim.plotAllStates(12, setdiff([1:24],[1, 15, 20, 11, 12, 7, 16, 21, 2]) );
% sim.plotAllOutputs(13,[1 32, 33, 34])


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
    



