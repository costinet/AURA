function inspect(obj)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    if ~isempty(obj.missingParams) 
        type = obj.missingParams.type;
        component = obj.missingParams.component;
        expression = obj.missingParams.val;

        T = cell2table([component', type', expression'], 'VariableNames', {'Component', 'Type', 'Expression'});
        disp(T);
        error(['Undefined parameters encountered when parsing PLECS circuit.  See table above for missing parameters. ', newline, ...
            'These parameters must be defined in the base workspace or altered in the PLECS circuit'])
    else
        try
            names = plecs('get', obj.sourcefn, 'Topology');
        catch e
            if startsWith(e.message, 'Error evaluating parameter')
                obj.readOnLoadError(e.message, obj.sourcefn);
                obj.inspect();
                return
            end
        end
    end

%     ssOrder = plecs('get', obj.sourcefn, 'StateSpaceOrder');
% 
%     modelFile = split(obj.sourcefn(), '/');
%     modelFile = join(modelFile(1:end-1), '/');
% 
%     allComponents = [ssOrder.States; ssOrder.Inputs; ssOrder.Switches];
%     for i = 1:length(allComponents)
%         componentPath = [modelFile{:} '/' allComponents{i}];
%         componentPath = split(componentPath,':');
%         componentPath = componentPath{1};
%         component = plecs('get', componentPath);
%         val = {};
%         switch component.Type
%             case 'Capacitor'
%                 val{1} = plecs('get', componentPath, 'C');
%             case 'Inductor'
%                 val{1} = plecs('get', componentPath, 'L');
%             case 'DCVoltageSource'
%                 val{1} = plecs('get', componentPath, 'V');
%             case 'Diode'
%                 val{1} = plecs('get', componentPath, 'Vf');
%                 val{2} = plecs('get', componentPath, 'Ron');
%             case 'MosfetWithDiode'
%                 val{1} = plecs('get', componentPath, 'Ron');
%             otherwise
%                 error(['Unknown component type ' component.Type]); 
%         end
% 
%         for j = 1:length(val)
%             try 
%                 numVal = evalin('base', val{j});
%             catch e
%                 disp(['Undefined parameter value "' val '" for ' component.Type ' ' component.Name ])
%             end
%         end
%     end
% 
%     try
%         names = plecs('get', obj.sourcefn, 'Topology');
%     catch e
%         
%     end
end