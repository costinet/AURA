function [T] = describeSwitchState(obj)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

    T = table;

    index = [1:length(obj.converter.swseq)]';
    swref = obj.converter.swseq';

    T = addvars(T, index, swref);

    for i = 1:size(obj.converter.swvec,2)
        T = addvars(T,  obj.converter.swvec(:,i));
        T.Properties.VariableNames{end} = obj.switchNames{i};
%         t(end) = obj.converter.swvec(:,i);
    end
    
    if nargout == 0
        disp(T)
    end
end