function linearizeCircuitModel(obj)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

    components = obj.origComponents;
   
    %% Replace FETs with Linear equivalent circuit
    FETs = strcmp({components.Type} , 'M');
    diodes = strcmp({components.Type} , 'D');
    switches = FETs | diodes;
    obj.topology.switchLabels = {obj.origComponents(switches).Name};

    numSwitches = sum(switches);
    numDiodes = sum(diodes);

    %% Replace Switches with linear subsystem
    locs = fliplr(find(switches));
    for i = 1:numSwitches 
        loc = locs(i);
        if FETs(loc)
            subCircuit = obj.switchLinearSubcircuit(components(loc));
        else
            assert(diodes(loc),'error');
            subCircuit = obj.diodeLinearSubcircuit(components(loc));
        end
        components = [components(1:loc-1), subCircuit, components(loc+1:end)];
    end


    %% Replace coupled inductors with equivalent circuit
    if ~isempty(obj.netListDirectives)
        for i = 1:length(obj.netListDirectives) 
            if obj.netListDirectives(i).Type == 'K'
                components = XFdependentSourceSubcircuit(obj, obj.netListDirectives(i), components);
            end
        end
    end


    %% Check for degenerate components that can be removed
    [components] = findShorts(obj, components);
    

    %% Store linearized component list
    obj.components = components;
   
end