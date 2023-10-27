function linearizeCircuitModel2(obj)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

    components = obj.origComponents;
    
    %% Replace FETs with Linear equivalent circuit
    FETs = strcmp({components.Type} , 'M');
    diodes = strcmp({components.Type} , 'D');
    switches = FETs | diodes;
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


    %% Replace coupled inductors with dependent sources
    if ~isempty(obj.netListDirectives)
        for i = 1:length(obj.netListDirectives) 
            if obj.netListDirectives(i).Type == 'K'
                components = XFdependentSourceSubcircuit(obj, obj.netListDirectives(i), components);
            end
        end
    end

%     %% Add measurement sources to everything else 
%     % This can eventually be deleted, but is currently used to match old
%     % implementation for comparison
%     otherPassives = [components.Type] == 'L' || [components.Type] == 'C';
%     locs = fliplr(find(otherPassives));
%     for i = 1:sum(otherPassives) 
%         loc = locs(i);
%         
%         components = [components(1:loc-1), subCircuit, components(loc+1:end)];
%     end



    obj.components = components;
   
end