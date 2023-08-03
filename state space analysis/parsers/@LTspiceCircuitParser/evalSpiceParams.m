function evalSpiceParams(obj,param)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    param = strrep(param,'{','');
    param = strrep(param,'}','');
    evalin('base', [obj.spiceNumFormat(param) ';']);

end