function component = parseSpiceComponent(obj, str)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    
    
    [sI, eI] = regexp(str, '[a-zA-Z][\S]*');
    name = strtrim(str(sI(1):eI(1)));

    type = name(1);

    % Find embedded params
    params = obj.parseSpiceParamList(str);

    component = struct;
    component.Name = name;
    component.Type = type;


    switch type 
        case 'M'
            %Need ron, roff, Coss
            vals = [0 obj.defaultRoff 0];

            [nodes,~,~] = regexp(str, '(?<=[\s])[\w]*(?=[\s]*)', 'match', 'tokens');
            component.Nodes = nodes([1 3]); % only drain and source used.

            [~, eI] = regexp(str, '[\w]*[\s]*[\w]*[\s]*[\w]*[\s]*[\w]*[\s]*[\w]*[\s]*'); %name, then 4 nets
            data = strtrim(str(eI(1):end));
            if length(data) <= 1
                %Can't be anything -- Ignore
            elseif ~contains(data, '=')
                % probably a part number
                if any(contains(obj.netlistModels,data))
                    sI = regexp(obj.netlistModels,['.model[\s]+',data,'[\s]+']);
                    modelLoc = ~cellfun(@isempty, sI);
                    params = obj.parseSpiceParamList(obj.netlistModels{modelLoc},params);
                else
                    % Look in obj.netlistLibraries
                    model = findModelInSpiceLibraries(obj,data);
                    params = obj.parseSpiceParamList(model,params);
                end
                vals = [.1 obj.defaultRoff 1e-9];
            end

            %params
            for i = 1:length(params)
                eqLoc = strfind( params{i}, '=');
                paramName = strtrim(params{i}(1:eqLoc-1));
                t = transistor();
                [paramName, isValid] = t.isParamOf(paramName);
%                 assert(isValid, ['Unknown transistor parameter "' params{i} '"'])
                if ~isValid
                    continue
                end


                paramVal = strtrim(params{i}(eqLoc+1:end));
                paramVal = spiceNumFormat(obj,paramVal);
                paramVal = strrep(paramVal,'{','');
                paramVal = strrep(paramVal,'}','');
                paramVal = evalin('base', paramVal);
                
                if strcmp(paramName, 'Rds')
                    vals(1) = paramVal;
                elseif strcmp(paramName, 'Coss') || strcmp(paramName, 'Cds')
                    vals(3) = paramVal;
                end

            end
            
            component.paramNames = {'Rds', 'Roff', 'Coss'};
            component.paramVals = vals;
        case 'V'
            [nodes,~,~] = regexp(str, '(?<=[\s])[\w]*(?=[\s]*)', 'match', 'tokens');
            component.Nodes = nodes([1 2]);

            component.paramNames = {type};
            component.paramVals = parseTwoNetSpiceComponent(obj,str);
        case 'I'
            [nodes,~,~] = regexp(str, '(?<=[\s])[\w]*(?=[\s]*)', 'match', 'tokens');
            component.Nodes = nodes([1 2]);

            component.paramNames = {type};
            component.paramVals = parseTwoNetSpiceComponent(obj,str);
        case 'L'
            [nodes,~,~] = regexp(str, '(?<=[\s])[\w]*(?=[\s]*)', 'match', 'tokens');
            component.Nodes = nodes([1 2]);

            component.paramNames = {type};
            component.paramVals = parseTwoNetSpiceComponent(obj,str);
        case 'C'
            [nodes,~,~] = regexp(str, '(?<=[\s])[\w]*(?=[\s]*)', 'match', 'tokens');
            component.Nodes = nodes([1 2]);

            component.paramNames = {type};
            component.paramVals = parseTwoNetSpiceComponent(obj,str);
        case 'R'
            [nodes,~,~] = regexp(str, '(?<=[\s])[\w]*(?=[\s]*)', 'match', 'tokens');
            component.Nodes = nodes([1 2]);

            component.paramNames = {type};
            component.paramVals = parseTwoNetSpiceComponent(obj,str);
        case 'D'
             %Need ron, roff, Cd, Vf
            vals = [0 obj.defaultRoff 0 0];

            [nodes,~,~] = regexp(str, '(?<=[\s])[\w]*(?=[\s]*)', 'match', 'tokens');
            component.Nodes = nodes([1 2]); 

            [~, eI] = regexp(str, '[\w]*[\s]*[\w]*[\s]*[\w]*[\s]*'); %name, then 2 nets
            data = strtrim(str(eI(1):end));
            if length(data) <= 1
                %Can't be anything -- Ignore
            elseif ~contains(data, '=')
                % probably a part number
                if any(contains(obj.netlistModels,data))
                    sI = regexp(obj.netlistModels,['.model[\s]+',data,'[\s]+']);
                    modelLoc = ~cellfun(@isempty, sI);
                    params = obj.parseSpiceParamList(obj.netlistModels{modelLoc},params);
                else
                    % Look in obj.netlistLibraries
                    model = findModelInSpiceLibraries(obj,data);
                    params = obj.parseSpiceParamList(model,params);
                end
                vals = [.1 obj.defaultRoff 1e-9];
            end

            %params
            for i = 1:length(params)
                eqLoc = strfind( params{i}, '=');
                paramName = strtrim(params{i}(1:eqLoc-1));
%                 t = transistor();
%                 [paramName, isValid] = t.isParamOf(paramName);
% %                 assert(isValid, ['Unknown transistor parameter "' params{i} '"'])
%                 if ~isValid
%                     continue
%                 end


                paramVal = strtrim(params{i}(eqLoc+1:end));
                paramVal = spiceNumFormat(obj,paramVal);
                paramVal = strrep(paramVal,'{','');
                paramVal = strrep(paramVal,'}','');
                paramVal = evalin('base', paramVal);
                
                if strcmpi(paramName, 'Vf')
                    vals(4) = paramVal;
                elseif ~strcmpi(paramName, 'Roff') && strcmpi(paramName(1), 'R')
                    vals(1) = paramVal;
                elseif strcmpi(paramName(1), 'C')
                    vals(3) = paramVal;
                end

            end
            
            component.paramNames = {'Ron', 'Roff', 'Cd', 'Vf'};
            component.paramVals = vals;
        case 'K'
            warning('Not sure how Jared wants this (K)')
        otherwise
            error(['Unrecognized netlist line ' str] );
    end

end