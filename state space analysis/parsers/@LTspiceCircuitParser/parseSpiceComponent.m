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
    component.paramExpressions = {};


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
                if ~isempty(obj.netlistModels) && any(contains(obj.netlistModels,data))
                    sI = regexp(obj.netlistModels,['.model[\s]+',data,'[\s]+']);
                    modelLoc = ~cellfun(@isempty, sI);
                    params = obj.parseSpiceParamList(obj.netlistModels{modelLoc},params);
                else
                    % Look in obj.netlistLibraries
                    model = findModelInSpiceLibraries(obj,data);
                    params = obj.parseSpiceParamList(model,params);
                end
                % default values if none other specified
                vals = [.1 obj.defaultRoff 1e-9];
                expr = {.1, obj.defaultRoff, 1e-9 };
            end

            %params
            for i = 1:length(params)
                %If any of the above generated lists of parameters, parse
                %them here to update values.
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
                try
                    paramExpr = paramVal;
                    paramVal = evalin('base', paramVal);
                catch e
%                      if (strcmp(e.identifier,'MATLAB:UndefinedFunction'))
%                       msg = ['Paramter ' paramName ' of component ' component.Name ' is not defined. ' ...
%                           'The expression ' paramVal ' is not defined in the netlist of in the base workspace.'];
%                         causeException = MException('MATLAB:LTSpiceParser:UndefinedParameter',msg);
%                         e = addCause(e,causeException);
%                     end
%                     rethrow(e)
                    obj.undefinedExpressions = [obj.undefinedExpressions; ...
                        {component.Name, paramName, paramExpr }];
                    paramVal = 0;
                    
                end
                
                if strcmp(paramName, 'Rds')
                    expr{1} = paramExpr;
                    vals(1) = paramVal;
                elseif strcmp(paramName, 'Coss') || strcmp(paramName, 'Cds')
                    expr{3} = paramExpr;
                    vals(3) = paramVal;
                end

            end
            
            component.paramNames = {'Rds', 'Roff', 'Coss'};
            component.paramVals = vals;
            component.paramExpressions = expr;
        case 'V'
            [nodes,~,~] = regexp(str, '(?<=[\s])[\w]*(?=[\s]*)', 'match', 'tokens');
            component.Nodes = nodes([1 2]);

            component.paramNames = {type};
            [component.paramVals, embeddedParams, component.paramExpressions] = parseTwoNetSpiceComponent(obj,str,component);
            component = deEmbedSpiceParams(obj,component,embeddedParams);
        case 'I'
            [nodes,~,~] = regexp(str, '(?<=[\s])[\w]*(?=[\s]*)', 'match', 'tokens');
            component.Nodes = nodes([1 2]);

            component.paramNames = {type};
            [component.paramVals, embeddedParams, component.paramExpressions] = parseTwoNetSpiceComponent(obj,str,component);
            component = deEmbedSpiceParams(obj,component,embeddedParams);
        case 'L'
            [nodes,~,~] = regexp(str, '(?<=[\s])[\w]*(?=[\s]*)', 'match', 'tokens');
            component.Nodes = nodes([1 2]);

            component.paramNames = {type};
            [component.paramVals, embeddedParams, component.paramExpressions] = parseTwoNetSpiceComponent(obj,str,component);
            component = deEmbedSpiceParams(obj,component,embeddedParams);
        case 'C'
            [nodes,~,~] = regexp(str, '(?<=[\s])[\w]*(?=[\s]*)', 'match', 'tokens');
            component.Nodes = nodes([1 2]);

            component.paramNames = {type};
            [component.paramVals, embeddedParams, component.paramExpressions] = parseTwoNetSpiceComponent(obj,str,component);
            component = deEmbedSpiceParams(obj,component,embeddedParams);
        case 'R'
            [nodes,~,~] = regexp(str, '(?<=[\s])[\w]*(?=[\s]*)', 'match', 'tokens');
            component.Nodes = nodes([1 2]);

            component.paramNames = {type};
            [component.paramVals, embeddedParams, component.paramExpressions] = parseTwoNetSpiceComponent(obj,str,component);
            component = deEmbedSpiceParams(obj,component,embeddedParams);
        case 'D'
             %Need ron, roff, Cd, Vf
            vals = [0 obj.defaultRoff 0 0];
            expr = {0, obj.defaultRoff, 0, 0};

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
                paramExpr = paramVal;
                paramVal = evalin('base', paramVal);
                
                if strcmpi(paramName, 'Vf')
                    vals(4) = paramVal;
                    expr{4} = paramExpr;
                elseif ~strcmpi(paramName, 'Roff') && strcmpi(paramName(1), 'R')
                    vals(1) = paramVal;
                    expr{41} = paramExpr;
                elseif strcmpi(paramName(1), 'C')
                    vals(3) = paramVal;
                    expr{3} = paramExpr;
                end

            end
            
            component.paramNames = {'Ron', 'Roff', 'Cd', 'Vf'};
            component.paramVals = vals;
            component.paramExpressions = expr;
        case 'K'
            [inductors,~,~] = regexp(str, '(?<=[\s])[\w]*(?=[\s]*)', 'match', 'tokens');

            component.paramNames = inductors(1:end-1);
            component.paramVals = inductors{end}; %actual K value
            component.paramExpressions = {};
            component.Nodes = {};
        otherwise
            warning(['Unrecognized netlist line ' str '.  This line will be ignored.'] );
    end

    for i = 1:length(component) %For loop in case embedded params de-embeeded
        if strcmp(component(i).Type,'V')
            if component(i).paramVals == 0
                % convert 0V sources to current measurement
                component(i).Type = 'Vm';
                component(i).Name(1) = 'I';
            end
        elseif strcmp(component(i).Type,'I')
            if component(i).paramVals == 0
                % convert 0A sources to voltage measurement
                component(i).Type = 'Im';
                component(i).Name(1) = 'V';
            end
        end
    end

end