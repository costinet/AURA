function [tLocs,insertAt,adjType] = findRequiredUncontrolledSwitching(obj,violateMarginStart,violateMarginEnd)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    [~,ints,~] = getIntervalts(obj.converter);
    ints = ints';

    %% "before" and "after" refer to time intervals
    errBefore = min(violateMarginStart,0)<0;  %There is an error at the start of the interval
    errAfter = min(violateMarginEnd,0)<0;     %There is an error at the end of the interval

    % Errors which start in one interval, but continue into the next
    interfaceErrors =  errBefore & circshift(errAfter,1,2);

    % Errors persisting for an entire interval
    fullIntErrors = errBefore & errAfter;

    % Errors persisting for an entire interval and into the next
    extendInterfaceErrors = fullIntErrors & circshift(errBefore,-1,2);

    endOfextendedError = extendInterfaceErrors & ~circshift(errAfter,-1,2);

    errLocs = (fullIntErrors & ~extendInterfaceErrors) | endOfextendedError ;
%     unadjustableLocs = (errAfter  & (ints ~= circshift(ints,-1))) | ...
%              (errBefore & (ints ~= circshift(ints,1)) & ~(circshift(errAfter,1,2)<0)) ;
    unadjustableLocs = (errAfter  & (ints ~= circshift(ints,-1))) | ...
            (errBefore & (ints ~= circshift(ints,1))) ;

    unadjustableLocs = unadjustableLocs & ~interfaceErrors;
    %These times cannot be adjusted because they are controlled switching
    %intervals.  The outer loop will need to insert a new uncontrolled
    %switching action before they can be altered.


    adjustableButSameInterval = (errAfter  & ...
        ((ints == circshift(ints,-1)) & (obj.converter.swind == circshift(obj.converter.swind,-1)))) | ...
        (errBefore & ...
        ((ints == circshift(ints,1)) & (obj.converter.swind == circshift(obj.converter.swind,1))));


    %% Three reasons to add a new interval
        % (1) there is error both before and after a single interface
        % (2) there is error either before or after, 
        %   (2a) but the time of the  interface isn't adjustable because 
        %         it is a controlled switching action
        %   (2b) but the before and after subinterval are the same (w.r.t.
        %   the error)
    tLocs = errLocs | unadjustableLocs | adjustableButSameInterval;
    

    %% Handle case where there are multiple connected intervals with violation 
    % because the violating switch doesn't change (i.e. others are switching)
%     tLocs(tLocs & circshift(tLocs,-1,2)) = 0; %If the following interval has the same error, don't also insert here
%             %start from the end and work your way back, instead
    

    insertAt = any(tLocs,1);
    
    adjType = zeros(size(obj.topology.Cbnd,1), length(obj.converter.ts), 2);
    adjType(:,:,1) = repmat(insertAt,size(errAfter,1),1) & (errAfter)+0;
    adjType(:,:,2) = ((-1)*(errBefore) & ~errAfter) & repmat(insertAt,size(errAfter,1),1);
end

