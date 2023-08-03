function str = spiceNumFormat(obj,str)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    if nargin == 0
        str = '1u';
    end
    
    try
        % If there aren't any non-numeric characters, just return the value
        val = eval(str);
        return
    catch
    end
    
    
    spiceVars = {'T' 'G', 'Meg', 'M', 'k', 'm', 'u', char(181), 'Âµ', 'n', 'p', 'f'};
    spiceVals = {'1e12' '1e9' '1e6' '1e6' '1e3' '1e-3' '1e-6' '1e-6' '1e-6' '1e-9' '1e-12' '1e-15'};

    matches = regexp(str,spiceVars);
    matches = ~cellfun(@isempty,matches);
    
    for i = find(matches)
        [sI,eI] = regexp(str,['(?<=[\-0-9.]*[0-9.])' spiceVars{i} '\>']);
        if ~isempty(sI)
            for j=length(sI):-1:1
                str =[str(1:sI(j)-1) '*' spiceVals{i} str(eI(j)+1:end)];
            end
        end
    
    end

end