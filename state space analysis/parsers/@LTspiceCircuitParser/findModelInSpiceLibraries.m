function model = findModelInSpiceLibraries(obj,data)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

    model = [];
    for i = 1:length(obj.netlistLibraries)
        fid=fopen(obj.netlistLibraries{i});
        text = char(fread(fid)');
        fclose(fid);
        lines = splitlines(text);

        sI = regexp(lines,['.model[\s]+',data,'[\s]+']);

        modelLoc = find(~cellfun(@isempty, sI));

        if ~isempty(modelLoc)
            model = lines{modelLoc};
            
%             % Find embedded params
%             [sI, eI] = regexp(line,'\s[a-zA-Z_]+[\s]*[=][\s]*[\w{}]+');
%             params = {};
%             if ~isempty(sI)
%                 for j=length(sI):-1:1
%                     params = [params; strtrim(str(sI(j):eI(j)))];
%                 end
%             end

        else
            continue
        end
    end

end

