function [paramVal, additionalParams] = parseTwoNetSpiceComponent(obj,str)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    [~, eI] = regexp(str, '[\w]*[\s]*[\w]*[\s]*[\w]*[\s]*'); %name, then 2 nets
    paramVal = strtrim(str(eI(1):end));
    paramVal = spiceNumFormat(obj,paramVal);
    paramVal = strrep(paramVal,'{','');
    paramVal = strrep(paramVal,'}','');

    % Look for any additional embedded params, e.g. "Rser=x"
    additionalParams = obj.parseSpiceParamList(paramVal);
    if ~isempty(additionalParams)
        for i = 1:length(additionalParams)
            paramVal = strrep(paramVal,additionalParams{i},'');
        end
    end

    paramVal = evalin('base', paramVal);

end