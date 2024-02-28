function [params] = parseSpiceParamList(obj, str, params)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
    
    if nargin == 2 
        params = {};
    end

    if isempty(str)
        return
    end
    
    [sI, eI] = regexp(str,'\s[a-zA-Z_]+[\s]*[=][\s]*[\w.{}()-*+/]+');
    if ~isempty(sI)
        for j=length(sI):-1:1
            params = [params; strtrim(str(sI(j):eI(j)))];
        end
    end
end