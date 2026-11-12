function [ts,swvec] = truncateMod(tsi,swveci, Tend)
%truncateMod Truncates a modulation waveform at a specified end time.
%
%   [ts, swvec] = truncateMod(tsi, swveci, Tend) truncates a sequence of 
%   time intervals (tsi) and their corresponding switching vector (swveci) 
%   so that the total time duration exactly matches Tend. 
%
%   Inputs:
%       tsi    - (1,N) numeric row vector. Represents the sequence of initial
%                time steps or durations for the modulation waveform.
%       swveci - (N,M) numeric matrix. Represents the switching states  
%                during each time step. The number of rows  
%                must exactly match the number of elements in tsi.
%       Tend   - (1,1) numeric scalar. The target end time at which to 
%                truncate the waveform. Must satisfy 0 <= Tend <= Ts
%
%   Outputs:
%       ts     - (1,K) numeric row vector. The truncated time steps. The sum 
%                of the elements in ts will exactly equal Tend.
%       swvec  - (K,M) numeric matrix. The truncated switching vectors 
%                corresponding to the time steps in ts.
%
%   See Also dutyMod, combineMods, plotModWf

arguments (Input)
    tsi (1,:) {mustBeNumeric,mustBeReal}
    swveci {mustBeNumeric,mustBeReal,mustMatchts(tsi,swveci)}
    Tend (1,1) {mustBeNumeric,mustBeReal}
end

% arguments (Output)
%     outputArg1
%     outputArg2
% end

    assert(Tend > 0 && Tend <= sum(tsi), 'Tend must occur between 0 and Ts of the modulation waveform')
    cumts = cumsum(tsi);
    ind = find(cumts > Tend,1,'first');
    if ind ==1
        ts = Tend;
        swvec = swveci(1,:);
    % elseif ind == length(tsi)
    %     ts = [tsi(1:ind-1), Tend-sum(tsi(1:ind-1))];
    %     swvec = swveci(1:ind,:);
    else
        if abs(Tend-sum(tsi(1:ind-1))) > 0
            ts = [tsi(1:ind-1), Tend-sum(tsi(1:ind-1))];
            swvec = swveci(1:ind,:);
        else
            ts = tsi(1:ind-1);
            swvec = swveci(1:ind-1,:);
        end
    end


end

function mustMatchts(ts,swvec)
if ~isequal(size(swvec,1),numel(ts))
    eid = 'Size:notEqual';
    msg = 'The number of rows of swvec must match the number of columns of ts';
    error(eid,msg)
end
end
