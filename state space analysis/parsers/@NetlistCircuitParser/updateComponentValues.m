function updateComponentValues(obj)
%updateComponentValues update all component values in the circuit from the
%base workspace
%
%   updateComponentValues(obj) for an LTSpiceCircuitParser object obj,
%   re-parses all state space matrices using currently available variables
%   from the base workspace where appropriate.  Compared to a complete
%   loadModel(obj) call, updateComponentValues will not re-read the netlist 
%   file, but solely update the values of components.  This saves time if
%   updating within a loop or optimizer.


    obj.undefinedExpressions = {};

    for i = 1:length(obj.origComponents)
        expr = obj.origComponents(i).paramExpressions;
        if ~iscell(expr)
            expr = {expr};
        end
    
        for j = 1:numel(expr)
            if ~isempty(expr{j})
                try
                    newParam = evalin('base',char(expr{j}));
                    if ~isempty(newParam)
                        obj.origComponents(i).paramVals(j) = newParam;
                    end
                catch e
                    obj.undefinedExpressions = [obj.undefinedExpressions; ...
                            {obj.origComponents(i).Name, '', char(expr{j}) }];
                    obj.origComponents(i).paramVals(j) = nan;
                end
            end
        end
    end

    obj.clearNumericalResults();

    if ~isempty(obj.undefinedExpressions)
        T = table(obj.undefinedExpressions(:,1), obj.undefinedExpressions(:,3), 'VariableNames',{'Component', 'Expression'});        
        warning(strjoin(["Parameter values undefined.  See table Below: ", newline,  formattedDisplayText(T)]));

        return;
    else
        obj.loadModel;

        % Try to populate input vector, if fully defined
        if isempty(obj.topology.converter.u) && ~isempty(obj.topology.labels.inputLabels)
            try 
                u= [];
                for i = 1:length(obj.topology.labels.inputLabels)
                    loc = strcmp({obj.origComponents.Name},obj.topology.labels.inputLabels(i));
                    assocComp = obj.origComponents(loc);
                    u(i) = assocComp.paramVals;
                end
                obj.topology.converter.u = u;
            catch
                % no action, just don't update u
            end
        end
    end
end