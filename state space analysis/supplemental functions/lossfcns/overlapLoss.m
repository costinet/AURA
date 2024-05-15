function [Pov, TurnOffCurrent, TurnOnVoltage]  = overlapLoss(sim,switchNames,switchModel, Ig)
%OVERLAPLOSS compute overlap loss from steady-state simulation
%   Detailed explanation goes here
    arguments
        sim SMPSim
        switchNames {mustBeText(switchNames)}
        switchModel {mustBeA(switchModel,'transistor'), ...
                mustMatchSize(switchNames,switchModel)}
        Ig {isnumeric, mustMatchSize(switchNames,Ig)}
    end

    if ~isa(switchNames,'cell')
        switchNames = {switchNames};
    end

    lastOn = zeros(1,size(sim.swvec,2));
    lastOff = zeros(1,size(sim.swvec,2));
    loopSwVec = [sim.swvec; sim.swvec(1,:)];
    [int,sw,~]= find(diff(loopSwVec,1) == -1);
    lastOn(sw) = int;
    [int,sw,~]= find(diff(loopSwVec,1) == 1);
    lastOff(sw) = int;

%     lastOn = mod(lastOn,size(sim.swvec,1));
%     lastOff = mod(lastOff,size(sim.swvec,1));


    
    Pov = zeros(numel(switchNames),1);
    TurnOffCurrent = zeros(numel(switchNames),1);
    TurnOnVoltage = zeros(numel(switchNames),1);
    for i = 1:numel(switchNames)
        sn = switchNames{i};
        swLoc = sim.sigLoc(sn,'sw');
        [Iname, Vname] = sim.parser.getSwitchMeasSourceNames(sn);
        idsLoc = sim.sigLoc(Iname,'y');
        vdsLoc = sim.sigLoc(Vname,'y');

        if lastOn(swLoc) == 0 || lastOff(swLoc) ==0
            Pov(i) = 0;
            continue
        end
        
        TurnOffCurrent(i) = sim.YsEnd(idsLoc,lastOn(swLoc));
        TurnOnVoltage(i) = sim.YsEnd(vdsLoc,lastOff(swLoc));

        if isscalar(switchModel)
            trans = switchModel;
        else
            trans = switchModel(i);
        end

        if isscalar(Ig)
            ig = Ig;
        else
            ig = Ig(i);
        end

        [~, ~, ~, ~, Qsw] = gateCharge(trans, 7);

        Pov(i) = Qsw/ig*TurnOnVoltage(i)*TurnOffCurrent(i)*1/2*1/sum(sim.ts);
        warning('This is currently incorrect, but is left for demonstration')
        % Need to know what the forced current is.
    end
end

function mustMatchSize(switchNames,switchModel)
% Used in argument validation
    msg = 'switchModel must either be one transistor, or one transistor for each switchName specified';
    if isa(switchNames,'cell')
        if ~isequal(numel(switchNames),numel(switchModel))  && ...
            numel(switchModel) ~= 1
            eid = 'Size:notEqual';            
            error(eid,msg)
        end
    else
        assert(numel(switchModel) == 1, msg);
    end
end