function [components] = deEmbedSpiceParams(obj,component,embeddedParams)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    components = component; 

    for i = 1:length(embeddedParams)
        eqLoc = strfind( embeddedParams{i}, '=');
        paramName = strtrim(embeddedParams{i}(1:eqLoc-1));
        paramVal = strtrim(embeddedParams{i}(eqLoc+1:end));
        paramVal = spiceNumFormat(obj,paramVal);
        paramVal = strrep(paramVal,'{','');
        paramVal = strrep(paramVal,'}','');
        expr = paramVal;
        try
            paramVal = evalin('base', paramVal);
        catch
            obj.undefinedExpressions = [obj.undefinedExpressions; ...
                {component.Name, paramName, expr }];
            paramVal = 0;
        end
                
        
        switch paramName
            case 'Rser'
                nodeToReplace = find(~strcmp(component.Nodes, '0'),1,'first');
                oldNode = component.Nodes{nodeToReplace};
                newNode = newIntNodeName(obj,component,nodeToReplace, 1);
                component.Nodes(nodeToReplace) = {newNode};

                newComponent = {};
                newComponent.Name = ['R_{' component.Name '_{Ser}}'];
                newComponent.Type = 'R';
                newComponent.Nodes = {oldNode, newNode};
                newComponent.paramExpressions = expr;
                newComponent.paramNames = {paramName};
                newComponent.paramVals = paramVal;

                components = [component, newComponent];
            case 'Rpar'
                warning(['Embedded parameter Rpar in component ' component.Name ' currently not supported in netlist embedded params; use a discrete component instead'])
            case 'Cpar'
                warning(['Embedded parameter Cpar in component ' component.Name ' currently not supported in netlist embedded params; use a discrete component instead'])
            case 'Lser'
                warning(['Embedded parameter Lser in component ' component.Name ' currently not supported in netlist embedded params; use a discrete component instead'])
            case 'RLshunt'
                warning(['Embedded parameter RLshunt in component ' component.Name ' currently not supported in netlist embedded params; use a discrete component instead'])
            case 'm'
                warning(['Embedded parameter m (device paralleling) in component ' component.Name ' currently not supported in netlist embedded params; explicitly alter component instead'])
            otherwise
                warning(['Unknown embedded parameter ' paramName ' in component ' component.Name ' will be ignored']);
        end
    end
end