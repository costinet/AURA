function [eta, Pout, Pin] = powerEff(obj,inComp,outComp)
%% CURRENTLY INCOMPLETE
arguments (Input)
    obj SMPSim
    inComp {isstring}
    outComp {isstring}
end

% arguments (Output)
%     outputArg1
%     outputArg2
% end
    Pin = zeros(numel(inComp),1);
    Pout = zeros(numel(outComp),1);

    if ~iscell(inComp)
        inComp = {inComp};
    end
    if ~iscell(outComp)
        inComp = {outComp};
    end

    for i = 1:max(numel(inComp), numel(outComp))
        if i <= numel(inComp)
            assert(any(strcmp(inComp{i}, {obj.parser.components.Name})), ...
                'inComp must be a cell array of component names that match components in the schematic');
            Pin(i) = findComponentPower(obj, inComp{i});
        end
        if i <= numel(outComp)
            assert(any(strcmp(outComp{i}, {obj.parser.components.Name})), ...
                'outComp must be a cell array of component names that match components in the schematic');
            Pout(i) = findComponentPower(obj, outComp{i});
        end
    end
end

function P = findComponentPower(obj, compName)
    comp = obj.parser.components(strcmp(compName, {obj.parser.components.Name}));
    if startsWith(compName, 'R')
    
    elseif startsWith(compName, 'V') 
        %prefer u over values
        if any(strcmp(compName, obj.inputNames))
            V = obj.u(strcmp(compName, obj.inputNames));
        else
            V = comp.paramVals;
        end
        %% BUT, GENERALLY WON'T HAVE CURRENT UNLESS MEASUREMENT ADDED
        C = obj.parser.getNewOutputVector(compName)
    
    elseif startsWith(compName, 'I')
    
    else
        error('Invalid component for input power')
    end
end