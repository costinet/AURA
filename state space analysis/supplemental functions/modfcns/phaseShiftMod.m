function [ts,swvec] = phaseShiftMod(ts,swvec, phase, propName, propVal)
%PHASESHIFTMOD shift modulation signals
%
%   [ts,swvec] = phaseShiftMod(ts,swvec, phase) uses phase in radians
%
%   [ts,swvec] = phaseShiftMod(ts,swvec, phase, 'propertyNAme', propertyVal)
%       properties include 'phaseUnits' which can be 'rad' (default),
%       'deg', or 'time'

    arguments
        ts (1,:) {mustBeNumeric,mustBeReal}
        swvec {mustBeNumeric,mustBeReal,mustMatchts(ts,swvec)}
        phase {mustBeNumeric,mustBeReal}
    end
    arguments (Repeating)
        propName {mustBeMember(propName,["phaseUnits"])}
        propVal {mustBeMember(propVal,["rad","deg", "time"])}
    end

    if isempty(propVal)
        phaseUnits = 'rad';
    else
        phaseUnits = propVal{strcmp(propName,'phaseUnits')};
    end

    Ts = sum(ts);

    if strcmp(phaseUnits, 'rad')
        tphi = mod(phase/(2*pi)*Ts,Ts);
    elseif strcmp(phaseUnits, 'deg')
        tphi = mod(phase/(360)*Ts,Ts);
    elseif strcmp(phaseUnits, 'time')
        tphi = mod(phase,Ts);
    end

    cumts = [cumsum(ts)] + tphi;
    rshift = 0;

    if tphi > 0
        brk = find(cumts > Ts, 1);
        rshift = sum(cumts > Ts);

        if ~isempty(brk)
            cumts = [cumts(1:brk-1), Ts, cumts(brk:end) - Ts];
            swvec = [swvec(1:brk-1,:); swvec(brk,:); swvec(brk:end,:)];

            swvec = circshift(swvec,rshift,1);    
            ts = diff([0 circshift(cumts,rshift)]); 
        end
    end

    zeroints = ts <= 10*max(eps(ts));
    ts(zeroints) = [];
    swvec(zeroints,:) = [];
end

function mustMatchts(ts,swvec)
    if ~isequal(size(swvec,1),numel(ts))
        eid = 'Size:notEqual';
        msg = 'The number of rows of swvec must match the number of columns of ts';
        error(eid,msg)
    end
end