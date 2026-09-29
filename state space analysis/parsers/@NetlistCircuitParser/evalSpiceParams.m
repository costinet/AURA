function evalSpiceParams(obj,param)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    param = strrep(param,'{','');
    param = strrep(param,'}','');
    
    % Replace whitespace/commas that separate parameter assignments with a single comma
    param = regexprep(param, '[\s,]+(?=[a-zA-Z_]\w*\s*=)', ', ');
    
    fmt = obj.spiceNumFormat(param);
    try
        evalin('base', [fmt ';']);
    catch ME
        disp(['Failed on string: ', fmt]);
        rethrow(ME);
    end
end