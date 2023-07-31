function [paramVal] = parseTwoNetSpiceComponent(obj,str)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    [~, eI] = regexp(str, '[\w]*[\s]*[\w]*[\s]*[\w]*[\s]*'); %name, then 2 nets
    paramVal = strtrim(str(eI(1):end));
    paramVal = spiceNumFormat(obj,paramVal);
    paramVal = strrep(paramVal,'{','');
    paramVal = strrep(paramVal,'}','');
    paramVal = evalin('base', paramVal);

end