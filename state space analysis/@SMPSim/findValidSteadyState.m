function niter = findValidSteadyState(obj)
%FINDVALIDSTEADYSTATE find steady-state solution to converter accounting for state-dependent switching actions.
%
%   FINDVALIDSTEADYSTATE(obj) takes an SMPSim object describing a switched
%   mode power converter and solves its steady-state operating waveforms
%   while accounting for state-dependent switching actions (such as diode
%   switching).  
%
%   The algorithm will take up to obj.maxItns iterations before returning
%   terminating.  If, on termination, the niter >= obj.maxItns, the
%   algorithm was not successful.  Otherwise, the SMPSim object will have a
%   valid steady-state solution loaded on termination.
%
%   See also SMPSim, steadyState

    niter = 1;

%     obj.steadyState;
    conv = obj.converter;
    top = conv.topology;

    if isempty(top.constraints.Cbnd)
        warning('No constraints specified.  In this case findValidSteadyState is equivalent to steadyState');
        obj.steadyState;
        return
    end

    %% finalRun 
    % once everything seems to be error-free based on discrete time points,
    % goes through once more with eigenvalue-based spacing to make sure no 
    % inter-sample violations are occuring.  
    finalRun = 0;
    
    %% Symmetry check (Experimental)
    % Checks for converters with half-cycle symmetry.  If one is found,
    % steady-state will only need to examine the first half-cycle
    if obj.allowHalfCycleReduction
        [TF,lastInt,Ihc] = obj.converter.checkForSymmetry;
        if(TF)
            obj.IHC = Ihc;
            swvec = conv.swvec;
            ts = conv.ts;
            top.loadCircuit(top.circuitParser.sourcefn,swvec(1:lastInt,:),1);
            conv.setSwitchingPattern(swvec(1:lastInt,:), ts(1:lastInt))
        end
    end
    
    %%  TimeStepping Attempt (Experimental)
    % If set, the algorithm will first time-step through one period looking
    % for diode switching before entering the steady-state solution.
    % Generally should not be used, as it is slow and requires the initial
    % times to be close such that one time-stepping period is 'close' to
    % steady-state
    if obj.timeSteppingInit > 0
        [Xf,ts,swinds,us] = timeSteppingPeriod(obj);
        
        numPeriods = obj.timeSteppingInit;
        Xsims = zeros(size(Xf,1),numPeriods);
        for i = 1:numPeriods
            conv.setSwitchingPattern(swinds, ts);
            conv.u = us;
            Xsims(:,i) = Xf;
            [Xf,ts,swinds,us] = obj.timeSteppingPeriod(Xf, ts, swinds, us);
            if obj.debug == 1
                disp([Xsims(:,i) Xf]);
            end
        end
    
        [newts,newswinds] = obj.format1DtimingVector(ts,swinds); 
        conv.u = us;
        conv.setSwitchingPattern(newswinds, newts)
        clear newts;
%         X0 = obj.Xs(:,1);
    end



while(1)
    obj.steadyState;

    %% Discrete timepoint violation margin
    [violateMarginStart,violateMarginEnd,targetValStart,targetValEnd] = obj.checkDiscreteErr;
    errBefore = min(violateMarginStart,0);
    errAfter = min(violateMarginEnd,0);

    if obj.debug2 == 1
        obj.describeDiscreteErrors;
    end


    %% Insert additional states as necessary
    % We only need to insert a new state if 
    %   -- The interface has an error both before and after the switching
    %   -- OR it has one of the above, and the before and after switching
    %   positions aren't part of the modifiable uncontrolled times.
    
    
    [tLocs,insertAt,adjType] = obj.findRequiredUncontrolledSwitching(violateMarginStart,violateMarginEnd);

    altered = 0;
    allChanges = [];

    switchRef = top.switchRef; %Important to store this now, as it could get altered inside the loop, below

    for i = flip(find(insertAt))
        [~, dts] = conv.getDeltaT();
        dt = max(min([min(dts),min(conv.controlledts)/10 , conv.ts(i)/10]), 1*conv.timingThreshold);
        for j = 1:2
            if(any(adjType(:,i,j)))
                [alt, newSwInd] = ...
                    conv.addUncontrolledSwitching(i,(-1)^(j+1), ...
                    dt,switchRef(tLocs(:,i),1),~switchRef(tLocs(:,i),2));
                altered = altered | alt;

                if obj.debug2 == 1
%                     interval, beforeAfter, newSwInd, switches, newStates
                    if altered ~= 0
                        locRefVec = [find(tLocs(:,i))]';
                        for locRef = locRefVec
                            if ~isempty(newSwInd)
                                allChanges = [allChanges; i, (-1)^(j+1), newSwInd, switchRef(locRef,1), ~switchRef(locRef,2)];
                            end
                        end
                    end
                end
            end
        end
    end

    if obj.debug2 == 1
        obj.describeInsertedIntervals(allChanges)
    end



    obj.steadyState;
    if(obj.debug)
        if ~exist('H','var')
            H = figure(obj.debugFigNo);
            H.Name = "Debug";
        end
        obj.plotAllStates(H);
        hold(H.Children,'on')
    end

    %% If anything was altered
    % in the above, we insert new switching sequences wherever needed to 
    % make sure we can (possibly) adjust timings to get zero error.  If we
    % changed the switching sequence, just loop back around until we have a
    % valid sequence.  Once the sequence isn't modified, proceed to adjust
    % timings.
    if ~altered
        if ~any(any(errBefore) | any(errAfter)) && finalRun == 0
            % If no errors in discrete time state vectors, run one more
            % time with higher time resolution to make sure no errors occur
            % in between discrete samples.
            finalRun = 1;
            if obj.finalRunMethod
                obj.steadyState;
                [~,ts,swinds] = obj.timeSteppingPeriod();
                conv.setSwitchingPattern(swinds, ts);
%                 obj.steadyState;
                if obj.debug2 == 1
                    disp('**Attempting Final Run with timeSteppingPeriod');
                end
            else
                eigs2tis(conv);
%                 obj.steadyState;
                if obj.debug2 == 1
                    disp('**Attempting Final Run with eigs2tis');
                end
            end
            continue;
        elseif ~any(any(errBefore) | any(errAfter)) && finalRun == 1
            % If no errors and we've done the above, we're finished
            break;
        else
            % Otherwise keep looping
            if finalRun ==1
                % If we previously inserted intervals through one of the
                % final run methods, delete out any unecessary ones now
                % that we've passed with ~altered (we should have the ones
                % we need)
                obj.converter.eliminateRedundtantTimeIntervals;
                finalRun = 0;
                continue;  
            end
            finalRun = 0;
            
        end
        
        %% Correct based on discrete jacobian
        [JoutStart,JoutEnd] = obj.discreteJacobianConstraint;
        
        [i1,j1] = find(errBefore);
        [i2,j2] = find(errAfter);

        A = zeros(length(j1)+length(j2),length(conv.ts));
        b = zeros(length(j1)+length(j2),1);
        e = zeros(length(j1)+length(j2),1);
%         tindex = 1:length(conv.ts);

        for i = 1:length(j1)
            %deltas is a vector of how the error constraint is affected at
            %the switching interface, by perturbations to every switching
            %interval's duration
            deltas = squeeze(JoutStart(i1(i),j1(i),:));
            A(i,:) = deltas';
            b(i,1) = targetValStart(i1(i),j1(i));
            e(i,1) = abs(targetValStart(i1(i),j1(i)) - violateMarginStart(i1(i),j1(i)));
        end
        if(isempty(i)), i=0; end
        for j = 1:length(j2)
            deltas = squeeze(JoutEnd(i2(j),j2(j),:));
            A(i+j,:) = deltas';
            b(i+j,1) = targetValEnd(i2(j),j2(j));
            e(i+j,1) = abs(targetValEnd(i2(j),j2(j)) - violateMarginEnd(i2(j),j2(j)));
        end


        % A is now (# of errors) x (number of time intervals) and 
        % b is (# of errors) x 1 

        unChangeable = isnan(sum(A,1));

        %unChangeable are generally the columns correspondint to controlled
        %timing intervals without any state-dependent switching actions in
        %them.

       

        %% Add zero net perturbation to controlled time intervals
        scaleF = norm(A(:,~unChangeable))/numel(A);
        [~, timeInts, ~] = conv.getIntervalts;
        A = [A; zeros(max(timeInts), size(A,2))];
        b = [b; zeros(max(timeInts), 1)];
        e = [e; conv.timingThreshold*ones(max(timeInts), 1)];
        for i = 1:max(timeInts)
            if ~any(timeInts'==i & unChangeable)
                A(end-max(timeInts)+i,timeInts==i) = scaleF;
            end
        end
        
        A(:,unChangeable) = [];
        emptyRows = all(A==0,2);
        b(emptyRows) = [];
        e(emptyRows) = [];
        A(emptyRows,:) = [];

        orig_state = warning;
        warning('off','MATLAB:singularMatrix')
        warning('off','MATLAB:nearlySingularMatrix')
        warning('off','MATLAB:rankDeficientMatrix')
        %         [msgstr, msgid] = lastwarn;

        tsolve = zeros(size(conv.ts));
        tsolve(~unChangeable) = -(A\b);


        %% Check for competing constraints
        % Speeds up DAB_Rload, as an example, but may be unnecessary
        % overall
        try
            numTimeRows = numel(unique(timeInts));
            nontimeRows = ones(length(emptyRows)-numTimeRows,1);
            nontimeRows(length(emptyRows),1) = 0;
    
            numTimeRows = sum(~nontimeRows&~emptyRows);
    
            constrDir = sign(A(1:end-numTimeRows,:));
            errTimeLocs = timeInts([j1;j2]);
            errTimeLocs = errTimeLocs(nontimeRows&~emptyRows);
            timeLocs = timeInts(~unChangeable);
    
            targets = [constrDir errTimeLocs];
            if any(targets(:,1:end-1) ==1,1) & any(targets(:,1:end-1) ==-1,1) 
                % There are competing constraints
                [X,Y] =  meshgrid(timeLocs,errTimeLocs);
                localChanges = (X==Y);
                localChanges(end+1:end+numTimeRows,:) = 1;
    
                localDir = sign(A.*localChanges);
                localDir = localDir(1:end-numTimeRows,:);
                [~,IA,IC] = unique(errTimeLocs);
                limLocalDir = localDir(IA,:);
                for i=2:numel(IC)
                    if IC(i)==IC(i-1)
                        limLocalDir(IC(i),:) = limLocalDir(IC(i),:)  .* (localDir(i,:) == localDir(i-1,:) );
                    end
                end
                localDir = limLocalDir(limLocalDir~=0)';
    
                solvedDir = sign(tsolve(~unChangeable));
                if numel(localDir) == numel(solvedDir)
                    if all(localDir ~= solvedDir)
                        A(~localChanges) = 0;
                        tsolve(~unChangeable) = -(A\b);
                    end
                end
            end
        catch
            1; %Do nothing.  set breakpoint for debug
        end

        %% Backup for unsuccessful solve
        if any(isnan(tsolve))
            tsolve(~unChangeable) = -pinv(A)*b;
%             warning('Sometimes this goes awry')            
        end

        

        %% Check whether solution is valid for local approximation
        solveWorked = abs(b+A*tsolve(~unChangeable)') < abs(e);
        solveErr = any(~solveWorked(e ~= conv.timingThreshold));
        while solveErr
            %Simultaneously satisfying all constraints is not possible
            %(see twoPhaseBuck for an example)
            err = abs(b+A*tsolve(~unChangeable)') - abs(e);
            [~,I] = max(err);
            I = find(err > .9*err(I));

            if isempty(I)
                error('failed to solve iteration in timing');
                break
            end

            A(I,:) = [];
            b(I,:) = [];
            e(I,:) = [];

            tsolve = zeros(size(conv.ts));
            tsolve(~unChangeable) = -(A\b);

            if any(isnan(tsolve))
                tsolve(~unChangeable) = -pinv(A)*b;
%                 warning('Sometimes this goes awry')            
            end

            solveWorked = abs(b+A*tsolve(~unChangeable)') < abs(e);
            solveErr = any(~solveWorked(e ~= conv.timingThreshold));
            
        end

        warning(orig_state);
        
        

        oldts = conv.ts;
        tps = conv.validateTimePerturbations2(tsolve) ;%/ convSpeed;

        conv.adjustUncontrolledTiming(1:length(tps), tps);        


        if obj.debug2 == 1
            obj.describeAlteredTimes(oldts);
        end
        
        if numel(obj.converter.ts) ~= numel(oldts)
            %If we lose an interval, check to see if the adjacent ones
            %should now be combined
            obj.converter.eliminateRedundtantTimeIntervals;
        end    
        

        if exist('newts','var')
            if numel(newts) == numel(conv.ts)
                if max(abs(newts - conv.ts)) < 10*conv.timingThreshold
                    % If the code gets in here, it looks like we're
                    % oscillating, so try to break out
                        if obj.finalRunMethod
                            obj.steadyState;
                            [~,ts,swinds] = obj.timeSteppingPeriod();
                            conv.setSwitchingPattern(swinds, ts);
                            obj.steadyState;
                            if obj.debug2 == 1
                                disp('Attempting to break oscilation with timeSteppingPeriod');
                            end
                        else
                            eigs2tis(conv);
                            obj.steadyState;
                            if obj.debug2 == 1
                                disp('Attempting to break oscilationwith eigs2tis');
                            end
                        end
                end
            end
        end

        newts = conv.ts;

        
        if(obj.debug)
            obj.steadyState;
            obj.plotAllStates(H);
%             hold(H.Children,'on')
%             for i = 1:length(H.Children)
%                 c = parula;
%                 lineCo = reshape([H.Children(i).Children.Color]',[3,length(H.Children(i).Children)])';
%                 ind = all(lineCo == [0    0.4470    0.7410],2);
%                 H.Children(i).Children(ind).Color = c((niter-1)*30+1,:);
%             end
        end
        
        niter = niter+1;
        
        if ~obj.suppressIterationOutput
            disp([niter sum(errBefore + errAfter, 'all')]);
        end

        
        if(niter > obj.maxItns)
            warning(['unable to solve valid Steady State within ' num2str(obj.maxItns) ' iterations'])
            break;
        end
        
        if(~any(tsolve))
            error('timing not modified');
        end

        

    end
end
end