function [params] = parseSpiceParamList(obj, str, params)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
    
    if nargin == 2 
        params = {};
    end

    if isempty(str)
        return
    end
    
%     [sI, eI] = regexp(str,'\s[a-zA-Z_]+[\s]*[=][\s]*[\w.{}()-*+/]+');
    [sI, eI] = regexp(str,'(?<=[\W])[a-zA-Z][\w]*[\s]*[=][\s]*[\w.{}()\-*+/]+');
    if ~isempty(sI)
        for j=length(sI):-1:1
            paramstr = strtrim(str(sI(j):eI(j)));
            % The final element in a library will catch the closing
            % parentheses, but we don't want to ignore all parentheses, so
            % check 
            if j==length(sI) && endsWith(paramstr,')')
                if sum(paramstr == ')') > sum(paramstr == '(')
                    paramstr = paramstr(1:end-1);
                end
            end
            params = [params; paramstr];
        end
    end
end