clear all;

niter = 0;
debug = 0;
debug2 = 0;
finalRunMethod = 0; %1 to use timestepping, 0 to use eigs2tis
timeSteppingInit = 0;

if ~debug
%     w = warning ('off','all');
else
    w = warning ('on','all');
end

sdir = mfilename('fullpath');
sdir = sdir(1:find(sdir=='\',1,'last')-1);
addpath(sdir);

%% Load test circuit
% modelfile = 'AsyncBoost'; PLECsModel = 'Boost_Async';
% modelfile = 'MRBuck'; PLECsModel = 'MRBuck';
% modelfile = 'DSC4to1'; PLECsModel = 'HDSC';
% modelfile = 'DAB'; PLECsModel = 'DAB_oneCap';
% modelfile = 'DABfull'; PLECsModel = 'DAB_8Cap';
modelfile = 'DSC4to1Diodesr2021a'; PLECsModel = 'HDSC_withDiodes';
% modelfile = 'Flyback'; PLECsModel = 'PC_Flyback';
% modelfile = 'SCBuckHybridBuck'; PLECsModel = 'BuckBuck';

% find_system(modelfile,'SearchDepth',1, 'IncludeCommented', 'on')
open_system(modelfile,'loadonly');
circuitPath = [modelfile '/' PLECsModel];
set_param(circuitPath,'Commented','on');
simout = sim(modelfile,eps);

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

Xss = sim.steadyState;
if(debug)
    sim.plotAllStates(1);
end

% ssOrder = plecs('get', circuitPath, 'StateSpaceOrder');
% outputs = ssOrder.Outputs;

outputs = top.outputLabels;

%% finalRun 
% once everything seems to be error-free based on discrete time points,
% goes through once more with eigenvalue-based spacing to make sure no 
% inter-sample violations are occuring.  
finalRun = 0;

%% Symmetry check
% May be useful but not doing anything with it yet.  Can identify that DAB,
% etc. exhibit half-cycle symmetry
TF = conv.checkForSymmetry;

%%  TimeStepping Attempt
if timeSteppingInit == 1
    tic
    [Xf,ts,swinds] = timeSteppingPeriod(sim);
    
    numPeriods = 500;
    Xsims = zeros(size(Xf,1),numPeriods);
    for i = 1:numPeriods
        warning('You cant do it this way, need to put them in the correct fullswinds structure');
        conv.setSwitchingPattern(swinds, ts);
    %     Xss = sim.steadyState;
        Xsims(:,i) = Xf;
        [Xf,ts,swinds] = sim.timeSteppingPeriod(Xf, ts, swinds );
        disp([Xsims(:,i) Xf]);
    end

    ts = [[ts(1); zeros(length(ts)-2,1)], [ts(2:end)]];
    swinds = [[swinds(1); zeros(length(swinds)-2,1)], [swinds(2:end)']];

%     load('flyBackTimeStepping')
    conv.setSwitchingPattern(swinds, ts)
    toc
    X0 = sim.Xs(:,1);
    % save('flyBackTimeStepping', 'Xf', 'swinds', 'swvec', 'ts')
end

tic
while(1)
    Xss = sim.steadyState;
    
    %% Update constraints per the current switching vector
%     [Cbnd, Dbnd, hyst, switchRef] = top.getPLECSConstraintMatrices(circuitPath);
    Cbnd = top.Cbnd; Dbnd = top.Dbnd; hyst = top.bndHyst; switchRef = top.switchRef;
    Cbnd = Cbnd(:,:,conv.swind);
    Dbnd = Dbnd(:,:,conv.swind);

    if(0)
        %% ~~Continuous time violation margin
        [ xs, t] = sim.SS_WF_Reconstruct;
        for i = 1:length(xs)-1
            swstate = find(t(i) <= cumsum(ts),1,'first');
            violationMargin(:,i) = Cbnd(:,:,swstate)*xs(:,i) + Dbnd(:,:,swstate)*us' - hyst(:,1) + hyst(:,2);
        end

        figure; plot(t(1:end-1),violationMargin)
        ylim([-100 100])
        legend(sim.switchNames) %% <- problem here, now.
    end

    %% Discrete timepoint violation margin
% %     violateMarginStart = zeros(size(Cbnd,1) ,length(conv.swind));
% %     targetValStart = violateMarginStart;
% %     violateMarginEnd = zeros(size(Cbnd,1) ,length(conv.swind));
% %     targetValEnd = violateMarginEnd;
% %     for i = 1:length(conv.swind)
% %         violateMarginStart(:,i) = Cbnd(:,:,i)*Xss(:,i) + Dbnd(:,:,i)*us' - hyst(:,1) + hyst(:,2);
% %         targetValStart(:,i) = Cbnd(:,:,i)*Xss(:,i) + Dbnd(:,:,i)*us' - hyst(:,1);
% %         violateMarginEnd(:,i) = Cbnd(:,:,i)*Xss(:,i+1) + Dbnd(:,:,i)*us' - hyst(:,1) + hyst(:,2);
% %         targetValEnd(:,i) = Cbnd(:,:,i)*Xss(:,i+1) + Dbnd(:,:,i)*us' - hyst(:,1);
% %     end
    [violateMarginStart,violateMarginEnd,targetValStart,targetValEnd] = sim.checkDiscreteErr;
    errBefore = min(violateMarginStart,0);
    errAfter = min(violateMarginEnd,0);

    %% Check if we have the correct state to change times

    % if yes, proceed to jacobian.
    % if no, add states, then redo steady-state & error detection.
    % OR -- just do it anyway, let addUncontrolledSwitching report on whether
    % it actually did anything.

    %% Insert additional states as necessary
    % We only need to insert a new state if 
    %   -- The interface has an error both before and after the switching
    %   -- OR it has one of the above, and the before and after switching
    %   positions aren't part of the modifiable uncontrolled times.
    
    [~,ints,~] = getIntervalts(conv);
    ints = ints';
%     if size(ints,1) > size(ints,2)
%         ints = ints';
%     end

%     errLocs = (errAfter<0 & circshift(errBefore,-1,2)<0);
%     unadjustableLocs = (errAfter<0 | circshift(errBefore,-1,2)<0) & (ints ~= circshift(ints,-1));
%     tLocs = errLocs | unadjustableLocs;
%     insertAfter = any(tLocs,1);
    % The above only worked for always inserting AFTER
    % Note, this was to not inser again when the error was addressable by
    % the altready
    
%     errLocs = (errAfter<0 & (errBefore<0));
%     unadjustableLocs = (errAfter<0  & (ints ~= circshift(ints,-1))) | ...
%         (errBefore<0 & (ints ~= circshift(ints,1)));
%     adjustableButSameInterval = (errAfter<0  & ...
%         ((ints == circshift(ints,-1)) & (conv.swind == circshift(conv.swind,-1)))) | ...
%         (errBefore<0 & ...
%         ((ints == circshift(ints,1)) & (conv.swind == circshift(conv.swind,1))));
%     tLocs = errLocs | unadjustableLocs | adjustableButSameInterval;
%     insertAt = any(tLocs,1);
%     
%     adjType = zeros(size(Cbnd,1), length(conv.ts), 2);
%     adjType(:,:,1) = (errAfter<0)+0;
%     adjType(:,:,2) = (-1)*(errBefore<0);
    
    [tLocs,insertAt,adjType] = sim.findRequiredUncontrolledSwitching(violateMarginStart,violateMarginEnd);

    altered = 0;
    
    
%     warning("Need something so it doesn't do corrections both before and after");
    for i = flip(find(insertAt))
        % addUncontrolledSwitching(obj, interval, beforeAfter, initialTime, switches, newStates)
        dt = min( min(conv.controlledts)/20 , conv.ts(i)/20);
        for j = 1:2
            if(any(adjType(:,i,j)))
                altered = altered | conv.addUncontrolledSwitching(i,(-1)^(j+1),dt,switchRef(tLocs(:,i),1),~switchRef(tLocs(:,i),2));
            end
        end
    end

    Xss = sim.steadyState;
    if(debug)
        sim.plotAllStates(1);
    end
    

    
    if ~altered
        
        if ~any(any(errBefore) | any(errAfter)) && finalRun == 0
            finalRun = 1;
            if finalRunMethod
                Xss = sim.steadyState;
                [Xf,ts,swinds] = sim.timeSteppingPeriod();
                conv.setSwitchingPattern(swinds, ts);
                Xss = sim.steadyState;
            else
                eigs2tis(conv);
                Xss = sim.steadyState;
            end
            continue;
        elseif ~any(any(errBefore) | any(errAfter)) && finalRun == 1
            break;
        else
            finalRun = 0;
        end
        %% Correct based on discrete jacobian
        
        [Jt, J2t] = discreteJacobian(sim, 2);
        
        %% Check validity of jacobian function (test)
        if(debug2)
                [Jo] = discreteJacobian(sim, 1);
                func = @(x) getSSforJacobian(sim, x);
                [jac,err] = jacobianest(func,zeros(size(conv.ts)));
                jac = reshape(jac, size(Jo))*1e12;
                Jo
                jac./Jo
%                 Jt = -jac
        end
   
        Jout = zeros(size(Cbnd,1),size(Jt,2),size(Jt,3));
        JoutStart = Jout;
        JoutEnd = Jout;
        %J(state, at time, time changed)
        for i = 1:size(Jt,3)
            for j = 1:size(Jt,2)
                Jout(:,j,i) = Cbnd(:,:,j)*Jt(:,j,i);

                % Cbnd is during the interval, but Jt is about the states at the
                % interface.
                JoutStart(:,j,i) = Cbnd(:,:,j)*Jt(:,j,i);
                endInt = circshift(1:size(Jt,2),-1);
                JoutEnd(:,j,i) = Cbnd(:,:,j)*Jt(:,endInt(j),i);
            end
        end
        
%         [JoutStart2,JoutEnd2] = sim.discreteJacobianConstraint

        
        %% Works with only one violation:
        % J = squeeze(Jout(:,intV,intV-1));
        % Err = violateMarginEnd
        % dt = -J\Err(:,intV);

% % %         A = zeros(length(conv.ts),length(conv.ts));
% % %         b = zeros(length(conv.ts),1);
% % %         tindex = 1:length(conv.ts);
% % % 
% % %         rows2elim = [];
% % %         cols2elim = [];
% % % 
% % %         for i = 1:length(j1)
% % %             deltas = squeeze(JoutEnd(i1(i),j1(i),:));
% % %             A(i,:) = deltas';
% % %             b(i,1) = violateMarginEnd(i1(i),j1(i));
% % %         end
% % %         if(isempty(i)) i=0; end
% % %         for j = 1:length(j2)
% % %             deltas = squeeze(JoutStart(i2(j),j2(j),:));
% % %             A(i+j,:) = deltas';
% % %             b(i+j,1) = violateMarginStart(i2(j),j2(j));
% % %         end
% % %         if(isempty(j)) j=0; end
% % %         for k = 1:(length(conv.ts) - (length(j1)+length(j2)))
% % %             candidateSkips = setdiff(1:length(conv.ts), ...
% % %                 [mod(j2-2,length(conv.ts))+1; j1]);
% % % 
% % %             rows2elim = [rows2elim, i+j+k];  %Equations 
% % %             cols2elim = [cols2elim, candidateSkips(k)]; % Times
% % % 
% % % 
% % %         %     A(i+j+k,:) = -NaN;
% % %         %     A(:,candidateSkips(k)) = -NaN;
% % %         %     b(i+j+k,:) = NaN;
% % %         %     tindex(candidateSkips(k)) = NaN;
% % %         end
% % % 
% % %         % A(isnan(A)) = [];
% % %         % b(isnan(b)) = [];
% % %         % tindex(isnan(tindex)) = [];
% % % 
% % %         A(rows2elim,:) = [];
% % %         A(:,cols2elim) = [];
% % %         b(rows2elim) = [];
% % %         tindex(cols2elim) = [];
% % % 
% % %         % for i = 1:length(intV)
% % %         %     deltas = squeeze(Jout(stateV(i),intV(i),:));
% % %         %     A(i,:) = deltas';
% % %         %     b(i,1) = min(violateMarginEnd(stateV(i),intV(i)), violateMarginStart(stateV(i),intV(i)));
% % %         % end
% % %         % for j = [length(intV):length(ts)] - (length(intV)-1)
% % %         %     candidateSkips = setdiff(1:length(ts), mod(intV-2,length(ts))+1);
% % %         %     A(i+j,candidateSkips(j)) = 1;
% % %         %     b(i+j,1) = 0;
% % %         % end

% % % %         [i1,j1] = find(violateMarginEnd < 0);
% % % %         [i2,j2] = find(violateMarginStart < 0);
% % % %         intV = [j1; j2]; stateV = [i1; i2];
% % % % 
% % % %         A = zeros(length(j1)+length(j2),length(conv.ts));
% % % %         b = zeros(length(j1)+length(j2),1);
% % % %         tindex = 1:length(conv.ts);
% % % % 
% % % %         for i = 1:length(j1)
% % % %             deltas = squeeze(JoutEnd(i1(i),j1(i),:));
% % % %             A(i,:) = deltas';
% % % %             b(i,1) = targetValEnd(i1(i),j1(i));
% % % %         end
% % % %         if(isempty(i)), i=0; end
% % % %         for j = 1:length(j2)
% % % %             deltas = squeeze(JoutStart(i2(j),j2(j),:));
% % % %             A(i+j,:) = deltas';
% % % %             b(i+j,1) = targetValStart(i2(j),j2(j));
% % % %         end
% % % %         if(isempty(j)), j=0; end
        
        [i1,j1] = find(errBefore);
        [i2,j2] = find(errAfter);
        intV = [j1; j2]; stateV = [i1; i2];

        A = zeros(length(j1)+length(j2),length(conv.ts));
        b = zeros(length(j1)+length(j2),1);
        tindex = 1:length(conv.ts);

        for i = 1:length(j1)
            deltas = squeeze(JoutStart(i1(i),j1(i),:));
            A(i,:) = deltas';
            b(i,1) = targetValStart(i1(i),j1(i));
        end
        if(isempty(i)), i=0; end
        for j = 1:length(j2)
            deltas = squeeze(JoutEnd(i2(j),j2(j),:));
            A(i+j,:) = deltas';
            b(i+j,1) = targetValEnd(i2(j),j2(j));
        end
        if(isempty(j)), j=0; end

        unChangeable = isnan(sum(A,1));
       
        
        tsolve = zeros(size(conv.ts));

        %% Attempt: add zero net perturbation to time as a part of the
        %% equations -- MAY NOT  because some times are dropped.
        scaleF = norm(A(:,~unChangeable))/numel(A);
        [~, timeInts, ~] = conv.getIntervalts;
        A = [A; zeros(max(timeInts), size(A,2))];
        b = [b; zeros(max(timeInts), 1)];
        for i = 1:max(timeInts)
            if ~any(timeInts'==i & unChangeable)
                A(end-max(timeInts)+i,timeInts==i) = scaleF;
            end
        end
        
        A(:,unChangeable) = [];

        tsolve(~unChangeable) = -(A\b);
        
        if any(isnan(tsolve))
            tsolve(~unChangeable) = -pinv(A)*b;
        end
        
        [~,dtLims] = getDeltaT(sim.converter);
        
        tr = tsolve./dtLims;
        tr = max(tr(tr>1));
        
        if ~isempty(tr)
            tsolve = tsolve/tr;
        end
        
%         deltaTs = zeros(size(conv.ts));
%         deltaTs(tindex) = tsolve;
        % 
        % %% for MRbuck initial iterations:
        % -violateMarginEnd(2,3)/Jt(3,4,3)

        %%HDSC only allow deadtimes to be adjusted
%         A(:,1) = [];
%         A(:,2) = [];
%         deltaDTs = -A\b
%         deltaTs = [0 deltaDTs(1) 0  deltaDTs(2)];

        oldts = conv.ts;
        for i=length(tsolve):-1:1
           if tsolve(i) ~= 0
%                conv.adjustTiming(i, deltaTs(i));
               conv.adjustUncontrolledTiming(i, tsolve(i)/2);
           end
        end
% % %         conv.adjustUncontrolledTimingVector(1:length(tsolve), tsolve)
        newts = conv.ts;

%         warning('Can get stuck trying to make intervals longer when it cannot');


        % for i = 1:length(intV)
        %     int = intV(i);
        %     J = squeeze(Jout(:,int,int-1));
        %     Err = violateMarginEnd;  %% THIS ISN't RIGHT, DEPENDS ON WHERE IT OCCURS
        %     dt(i) = -J\Err(:,int);
        %     conv.adjustTiming(intV(i), dt(i));
        % end


        % % Cbnd(:,:,3)*Jt(:,4,3)
        % % dt = -  Cbnd(:,:,3)*Xss(:,4) [2]   /  Cbnd(:,:,3)*Jt(:,4,3) [2]
        % %     =   -Err  /  Jacobian (assuming a state to correct from)
        % sim.ts(3) = sim.ts(3) + dt;
        % sim.ts(4) = sim.ts(4) - dt;
        
        if(debug)
%             disp(tsolve)
%             disp(sum(errBefore + errAfter, 'all'))
            Xss = sim.steadyState;
            sim.plotAllStates(10);
        end
        
        niter = niter+1;
        disp([niter sum(errBefore + errAfter, 'all')]);
        
        if(~any(tsolve))
            error('timing not modified');
            break;
        end

    end
end
toc

Xss = sim.steadyState;
sim.plotAllStates(10);
% sim.plotAllStates(11, [1, 15, 20, 11, 12, 7, 16, 21, 2]  );
% sim.plotAllStates(12, setdiff([1:24],[1, 15, 20, 11, 12, 7, 16, 21, 2]) );
% sim.plotAllOutputs(13,[1 32, 33, 34])

function Xs = getSSforJacobian(sim, newTs)
    [tps] = sim.converter.validateTimePerturbations([1:length(newTs)], newTs/1e12);
    Xs = sim.perturbedSteadyState(tps);
    Xs = Xs(:,1:end-1);
end



