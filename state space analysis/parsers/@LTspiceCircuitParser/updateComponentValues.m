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
                newParam = evalin('base',char(expr{j}));
                if ~isempty(newParam)
                    obj.origComponents(i).paramVals(j) = newParam;
                end
            end
        end
    end

    obj.Anum = [];
    obj.components = [];

    obj.loadModel;
end