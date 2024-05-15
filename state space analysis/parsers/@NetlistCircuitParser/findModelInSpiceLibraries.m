function model = findModelInSpiceLibraries(obj,data)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

    model = [];
    for i = 1:length(obj.netlistLibraries)
        fid=fopen(obj.netlistLibraries{i});
        if fid == -1
            % check if this is a problem with host default LTspice location
            fnRem =  strfind(obj.netlistLibraries{i},'LTspiceXVII');
            if fnRem
                subDir = obj.netlistLibraries{i}(fnRem+12:end);
                fid=fopen([obj.LTSpiceFolder '\' subDir]);
            end

            if fid == -1
                warning(['Unable to locate Spice library ' obj.netlistLibraries{i} ...
                    '.  Edit the provided netlist or try setting LTSpiceFolder in @LTSpiceCircuitParser']);
                continue
            end
        end
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
    if isempty(model)
        warning(['Unable to find definition for model ' data '.  Will attempt to continue with default parameters for the component'])
    end

end

