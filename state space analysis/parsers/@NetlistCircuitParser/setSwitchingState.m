function setSwitchingState(obj, swvec)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

    %% Find all switches in original (pre-linearization) circuit
    FETs = strcmp({obj.origComponents.Type} , 'M');
    diodes = strcmp({obj.origComponents.Type} , 'D');
    switches = FETs | diodes;
    numSwitches = sum(switches);
    numDiodes = sum(diodes);

    
    locs = fliplr(find(switches));
    for i = 1:numSwitches 
        loc = locs(i);

        % For each switch, find the resistor it was replaced with and its 
        % position in the swvec
        switchResistorLoc = strcmp({obj.components.Name}, ['R_' obj.origComponents(loc).Name]);
        swvecLoc = strcmp(obj.topology.switchLabels, [obj.origComponents(loc).Name]);

        % replace resistor's resistance with on/off value as appropriate
        if swvec(swvecLoc) == 1
            obj.components(switchResistorLoc).paramVals = ...
                obj.origComponents(loc).paramVals( ...
                strcmp(obj.origComponents(loc).paramNames,'Rds') | ...
                strcmp(obj.origComponents(loc).paramNames,'Ron'));
        else
            obj.components(switchResistorLoc).paramVals = ...
                obj.origComponents(loc).paramVals( ...
                strcmp(obj.origComponents(loc).paramNames,'Roff'));
        end
       
    end
end