function [tLocs,insertAt,adjType] = findRequiredUncontrolledSwitching(obj,violateMarginStart,violateMarginEnd)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    [~,ints,~] = getIntervalts(obj.converter);
    ints = ints';

    errBefore = min(violateMarginStart,0);
    errAfter = min(violateMarginEnd,0);

    errLocs = (errAfter<0 & (errBefore<0));
    unadjustableLocs = (errAfter<0  & (ints ~= circshift(ints,-1))) | ...
        (errBefore<0 & (ints ~= circshift(ints,1)));
    adjustableButSameInterval = (errAfter<0  & ...
        ((ints == circshift(ints,-1)) & (obj.converter.swind == circshift(obj.converter.swind,-1)))) | ...
        (errBefore<0 & ...
        ((ints == circshift(ints,1)) & (obj.converter.swind == circshift(obj.converter.swind,1))));
    tLocs = errLocs | unadjustableLocs | adjustableButSameInterval;
    insertAt = any(tLocs,1);
    
    adjType = zeros(size(obj.topology.Cbnd,1), length(obj.converter.ts), 2);
    adjType(:,:,1) = (errAfter<0)+0;
    adjType(:,:,2) = (-1)*(errBefore<0);
end

