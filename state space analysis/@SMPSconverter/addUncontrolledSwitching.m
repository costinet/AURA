function [altered, newSwInd] = addUncontrolledSwitching(obj, interval, beforeAfter, initialTime, switches, newStates, force)
    
    try 
        assert(initialTime < obj.ts(interval), 'specified initialTime is larger than the interval it subdivides');
    catch
        initialTime = obj.ts(interval)/2;
    end
    
    if nargin < 7
        force = 0;
    end

    altered = 0;
    newSwInd  = [];



    if initialTime < .9*obj.timingThreshold
%         warning('should this be here?')
        obj.eliminateZeroTimeIntervals
        return
    end

%     [r,c] = find(obj.fullts ~=0);
%     if size(r,1) < size(r,2)
%         r = r';
%         c = c';
%     end
%     nonZeroLocs = [r,c];
% 
%     target = nonZeroLocs(interval,:);

%     if obj.fullts(end,target(2)) == 0
%         replace = 1:size(obj.fullts,1);
%     else
%         replace = 1:size(obj.fullts,1)+1;
%     end
% 
%     % Make sure topology has the new switch vector parsed
%     newSwVec = obj.swvec(interval,:);
%     newSwVec(switches) = newStates;
%     if ~ismember(newSwVec, obj.topology.swseq, 'rows')
%         %with force == 0, should just add the additional switching
%         %states
%         obj.topology.loadCircuit(obj.topology.sourcefn, [obj.swvec; newSwVec], 0);
%     end
%     [~, newSwInd] = ismember(newSwVec, obj.topology.swseq, 'rows');
%     
%     %% Check for duplicate states
%     inds = find(nonZeroLocs(:,2) == target(2));
%     targetLoc = find(inds == interval);
%     
%     
%     altered = 0;
%     if targetLoc == 1 && beforeAfter == -1
%         altered = ~(obj.swind(interval) == newSwInd);
%     elseif targetLoc == length(inds) && beforeAfter == 1
%         altered = ~(obj.swind(interval) == newSwInd);
%     elseif  beforeAfter == -1
%         altered = ~(obj.swind(interval) == newSwInd || ...
%             (obj.swind(interval-1) == newSwInd) && 1    );
%     elseif beforeAfter == 1
%         altered = ~(obj.swind(interval) == newSwInd || obj.swind(interval+1) == newSwInd);
%     else
%         error('invalid case');
%     end
%     
%     if ~altered
%         warning('No state insertion performed because state already exists in this location');
%         return
%     end
%     
% 
%     if beforeAfter == -1 %before
%         obj.fullts(replace,target(2)) = ...
%            [obj.fullts(1:target(1)-1,target(2)); 
%            initialTime;
%            obj.fullts(target(1),target(2)) - initialTime;
%            obj.fullts(target(1)+1:replace(end)-1,target(2))]; 
%         obj.swind = [obj.swind(1:interval-1), newSwInd, obj.swind(interval:end)];
%     elseif beforeAfter == 1%after
%         obj.fullts(replace,target(2)) = ...
%            [obj.fullts(1:target(1)-1,target(2)); 
%            obj.fullts(target(1),target(2)) - initialTime;
%            initialTime;
%            obj.fullts(target(1)+1:replace(end)-1,target(2))]; 
%         obj.swind = [obj.swind(1:interval), newSwInd, obj.swind(interval+1:end)];
%     else
%         error('input variable beforeAfter must be +1 or -1');
%     end
%     
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Find interval and subinterval of target location
    [ts, ints, subInts] = getIntervalts(obj);
    oint = ints(interval);                  
    subInt = subInts(interval);
    
    if obj.fullts(end,oint) == 0
        replace = 1:size(obj.fullts,1);
    else
        replace = 1:size(obj.fullts,1)+1;
    end

    %% Make sure topology has the new switch vector parsed
    newSwVec = obj.swvec(interval,:);
    newSwVec(switches) = newStates;
    if ~ismember(newSwVec, obj.topology.swseq, 'rows')
        %with force == 0, should just add the additional switching
        %states
%         obj.topology.loadCircuit(obj.topology.sourcefn, [obj.swvec; newSwVec], 0);
        obj.topology.loadCircuit([], [obj.swvec; newSwVec], 0);
    end
    [~, newSwInd] = ismember(newSwVec, obj.topology.swseq, 'rows'); 
    
    %% Check for duplicate states
    altered = 0;
    if interval == 1 && beforeAfter == -1                   %before the first interval
        altered = ~(obj.swind(interval) == newSwInd);       %... if it isn't already there
    elseif interval == length(obj.ts) && beforeAfter == 1  %after the last interval
        altered = ~(obj.swind(interval) == newSwInd);       %... if it isn't already there
    elseif  beforeAfter == -1                               %somewhere in the middle
        altered = ~(obj.swind(interval) == newSwInd || ...  %... if it isn't there and is alterable
            ((obj.swind(interval-1) == newSwInd) && (ints(interval-1) == oint)));
    elseif beforeAfter == 1                                 %somewhere in the middle
        altered = ~(obj.swind(interval) == newSwInd || ...  %... if it isn't there and is alterable
            ((obj.swind(interval+1) == newSwInd) && (ints(interval+1) == oint)));
    else
        error('invalid case');
    end
    
    if ~altered && ~force
%         warning('No state insertion performed because state already exists in this location');
        return
    end
    
    %% Otherwise, do the replacement
    if beforeAfter == -1 %before
        obj.fullts(replace,oint) = ...
           [obj.fullts(1:subInt-1,oint); 
           initialTime;
           obj.fullts(subInt,oint) - initialTime;
           obj.fullts(subInt+1:replace(end)-1,oint)]; 
%         obj.swind = [obj.swind(1:interval-1), newSwInd, obj.swind(interval:end)];
         obj.fullswind(replace,oint) = ...
           [obj.fullswind(1:subInt-1,oint); 
           newSwInd;
           obj.fullswind(subInt:replace(end)-1,oint)];
    elseif beforeAfter == 1%after
        obj.fullts(replace,oint) = ...
           [obj.fullts(1:subInt-1,oint); 
           obj.fullts(subInt,oint) - initialTime;
           initialTime;
           obj.fullts(subInt+1:replace(end)-1,oint)]; 
%         obj.swind = [obj.swind(1:interval), newSwInd, obj.swind(interval+1:end)];
        obj.fullswind(replace,oint) = ...
           [obj.fullswind(1:subInt,oint); 
           newSwInd;
           obj.fullswind(subInt+1:replace(end)-1,oint)];
    else
        error('input variable beforeAfter must be +1 or -1');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    altered = 1;
    assert(all(abs((sum(obj.fullts,1) - obj.controlledts))<obj.timingThreshold), ...
        'Timing mismatch due to uncontrolled switching');

    obj.limitTimeIntervals();
end
