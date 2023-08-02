function varargout = findValidSteadyState(obj)
%FINDVALIDSTEADYSTATE is function of class SMPSim and corrects steady-state
%errors caused by diodes
%   FINDVALIDSTEADYSTATE() is called to eliminate errors cause by improper
%   dependent switching. 
%
%   By setting debug2 = 1 the describe steady-state  

    debug = 0;
    debug2 = 0;

    timeSteppingInit = 0;
    finalRunMethod = 0;
    niter = 1;
    allAssess = 1;

    tmult = 0.5;

    Xss = obj.SS_Soln();
    conv = obj.converter;
    top = conv.topology;
    trLimit = 10;
    trTracker = 0;

    %% finalRun 
    % once everything seems to be error-free based on discrete time points,
    % goes through once more with eigenvalue-based spacing to make sure no 
    % inter-sample violations are occuring.  
    finalRun = 0;
    
    %% Symmetry check
    % May be useful but not doing anything with it yet.  Can identify that DAB,
    % etc. exhibit half-cycle symmetry
%     TF = obj.converter.checkForSymmetry;
    
    %%  TimeStepping Attempt
    if timeSteppingInit > 0
        [Xf,ts,swinds] = timeSteppingPeriod(obj);
        
        numPeriods = timeSteppingInit;
        Xsims = zeros(size(Xf,1),numPeriods);
        for i = 1:numPeriods
            conv.setSwitchingPattern(swinds, ts);
            Xsims(:,i) = Xf;
            [Xf,ts,swinds] = obj.timeSteppingPeriod(Xf, ts, swinds );
            if debug == 1
                disp([Xsims(:,i) Xf]);
            end
        end
    
        [newts,newswinds] = obj.format1DtimingVector(ts,swinds);
        
    
    
        conv.setSwitchingPattern(newswinds, newts)
        clear newts;
        X0 = obj.Xs(:,1);
    end

    

    

while(1)
    Xss = obj.SS_Soln();
    
    if debug2
        obj.converter.swvec
        obj.converter.ts
        %obj.plotAllStates(niter)
    end

    %% Update constraints per the current switching vector
    Cbnd = top.Cbnd; Dbnd = top.Dbnd; 
    hyst = top.bndHyst; switchRef = top.switchRef;
    Cbnd = Cbnd(:,:,conv.swind);
    Dbnd = Dbnd(:,:,conv.swind);

    %% Discrete timepoint violation margin
    [violateMarginStart,violateMarginEnd,targetValStart,targetValEnd] = obj.checkDiscreteErr;
    errBefore = min(violateMarginStart,0);
    errAfter = min(violateMarginEnd,0);




    if debug2 == 1
        obj.describeDiscreteErrors;
    end




    %% Insert additional states as necessary
    % We only need to insert a new state if 
    %   -- The interface has an error both before and after the switching
    %   -- OR it has one of the above, and the before and after switching
    %   positions aren't part of the modifiable uncontrolled times.
    
    [~,ints,~] = getIntervalts(conv);
    ints = ints';
    
    [tLocs,insertAt,adjType] = obj.findRequiredUncontrolledSwitching(violateMarginStart,violateMarginEnd);

    % Aded by Baxter to only look at the end of times initally as it is
    % where oftne there are the most errors to correct
    %%{
    if niter < 20
        if finalRun == 0
            adjType(:,:,2) = zeros(size(adjType(:,:,2)));
            finalRun = 0;
        end
    end
    %}
    altered = 0;
    allChanges = [];
%     warning("Need something so it doesn't do corrections both before and after");

    for i = flip(find(insertAt))
        % addUncontrolledSwitching(obj, interval, beforeAfter, initialTime, switches, newStates)
%         dt = min( min(conv.controlledts)/20 , conv.ts(i)/20);
        [~, dts] = conv.getDeltaT();
        dt = min([min(dts),min(conv.controlledts)/20 , conv.ts(i)/20]);
        for j = 1:2
            if(any(adjType(:,i,j)))
                [alt, newSwInd] = ...
                    conv.addUncontrolledSwitching(i,(-1)^(j+1), ...
                    dt,switchRef(tLocs(:,i),1),~switchRef(tLocs(:,i),2));
                altered = altered | alt;

                if debug2 == 1
%                     interval, beforeAfter, newSwInd, switches, newStates
                    locRefVec = [find(tLocs(:,i))]';
                    for locRef = locRefVec
                        allChanges = [allChanges; i, (-1)^(j+1), newSwInd, switchRef(locRef,1), ~switchRef(locRef,2)];
                    end
%                     if i == find(insertAt,1,'first')
%                         obj.describeInsertedIntervals(allChanges)
%                     end
                end
            end
        end
    end

    if debug2 == 1
        obj.describeInsertedIntervals(allChanges)
    end


    Xss = obj.SS_Soln;
    if(debug)
        if ~exist('H','var')
            H = figure;
            H.Name = "Debug";
        end
        obj.plotAllStates(H);
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
            if finalRunMethod
                Xss = obj.steadyState;
                [Xf,ts,swinds] = obj.timeSteppingPeriod();
                conv.setSwitchingPattern(swinds, ts);
                Xss = obj.steadyState;
            else
                eigs2tis(conv);
                Xss = obj.steadyState;
            end
            continue;
        elseif ~any(any(errBefore) | any(errAfter)) && finalRun == 1
            % If no errors and we've done the above, we're finished
            break;
        else
            % Otherwise keep looping
            finalRun = 0;
        end
        
        %% Correct based on discrete jacobian
        [JoutStart,JoutEnd] = obj.discreteJacobianConstraint;
        
        [i1,j1] = find(errBefore);
        [i2,j2] = find(errAfter);
        intV = [j1; j2]; 
        stateV = [i1; i2];

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

        %% Check for slope change (Nonfunctional)
%         interinterval = (errAfter~= 0 & errBefore == 0);
%         sign(JoutStart) ~= sign(JoutEnd);

       

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

        tsolve = zeros(size(conv.ts));
        tsolve(~unChangeable) = -(A\b);
        
        if any(isnan(tsolve))
            tsolve(~unChangeable) = -pinv(A)*b;
            warning('Sometimes this goes awry')

%             tsolve = tsolve*min(conv.ts)/max(abs(tsolve))*.01;
            
        end

        % BAXTER This a added to not change values that are right next to errors
        % that are found where there are states that are eliminated even if
        % there is not an error there
        unChangeable2  = unChangeable|~(circshift(sum(errAfter~=0),1)~=0  | sum(errAfter~=0)~=0  | circshift(sum(errAfter~=0),-1)~=0);

       % if unChangeable2

        oldts = conv.ts;


        [~,dtLims] = getDeltaT(obj.converter);

        tr = tsolve./dtLims;
        tr = max(tr(tr>1));

        if tr>10
            J = 4564564;
        end

        if ~isempty(tr)
            if trTracker <trLimit
            tsolve = tsolve/tr;
            trTracker =trTracker+1;
            else
                trTracker = 0;
            end

        else
            trTracker = 0;
        end

        tps = conv.validateTimePerturbations2(tsolve);
        for i=length(tsolve):-1:1
            if tsolve(i) ~= 0
                conv.adjustUncontrolledTiming(i, tps(i));
            end
        end

        if debug2 == 1
            obj.describeAlteredTimes(oldts);
        end
%         
%         [~,dtLims] = getDeltaT(obj.converter);
%         
%         tr = tsolve./dtLims;
%         tr = max(tr(tr>1));
%         
%         if ~isempty(tr)
%             tsolve = tsolve/tr;
%         end
%         
        

        
%         recurs = 0;
%         for i=length(tsolve):-1:1
%            if tsolve(i) ~= 0
% %                conv.adjustTiming(i, deltaTs(i));
%                 
%                 %% If adjustUncontrolledTiming has an issue and changes the interval
%                 % should use validateTimePerturbations to find this on the
%                 % front end
%                 if recurs ~=1
%                     [~, recurs] = conv.adjustUncontrolledTiming(i, tsolve(i));
%                 else
%                     recurs = 0;
%                 end
%            end
%         end
% % %         conv.adjustUncontrolledTimingVector(1:length(tsolve), tsolve)
%{       
if exist('newts','var')
            if numel(newts) == numel(conv.ts)
                if max(abs(newts - conv.ts)) < 10*conv.timingThreshold
                    % If the code gets in here, it looks like we're
                    % oscillating, so try to break out
                    [Xf,ts,swinds] = timeSteppingPeriod(obj);
                    [newts,newswinds] = obj.format1DtimingVector(ts,swinds);
                    conv.setSwitchingPattern(newswinds, newts)
                    warning('hack')
                end
            end
        end
%}
        newts = conv.ts;

%         [[newts-oldts]', tsolve']



%         warning('Can get stuck trying to make intervals longer when it cannot');

        
        if(debug)
%             disp(tsolve)
%             disp(sum(errBefore + errAfter, 'all'))
            Xss = obj.steadyState;
            obj.plotAllStates(H);
        end
        
        niter = niter+1;
        
%         if(~allAssess)
% %             weightTotalErr = getWeightedTotalError(obj, errBefore,errAfter);
% %             disp([niter weightTotalErr]);
            disp([niter sum(errBefore + errAfter, 'all')]);
%         else
% %             updateWaitBar(h, modelfile, selectedModel, niter, altered, finalRun);
%         end
        
        if(allAssess && niter > 250)
            break;
        end
        
        if(~any(tsolve))
            error('timing not modified');
            break;
        end

        

    end
end

if(debug)
    if ~exist('H','var')
        close(H)
    end
end

varargout{1} = niter;
nargout = 1;
end