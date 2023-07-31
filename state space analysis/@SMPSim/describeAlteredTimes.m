function describeAlteredTimes(obj,oldts)
%describeAlteredTimes Print to console a description of how timing
%intervals were altered
%   
%   describeAlteredTimes(obj,oldts) 
%   prints out a description of how an SMPSim objects current timing intervals 
%   differ from oldts.  Results are only calculated if the same number of 
%   subintervals are present in both the current and old time vector.
%   Used for debugging only.
%
%   See Also SMPSim.describeDiscreteErrors,
%   SMPSim.describeInsertedIntervals, SMPSim.describeSwitchState

if all(size(obj.ts) == size(oldts))

    deltaTs = obj.ts - oldts;
    
    if all(deltaTs == 0)
        disp('No changes made to interval timings');
    else
        for i = 1:length(deltaTs)
            if deltaTs(i) < 0
                disp(['Duration of subinterval ' num2str(i) ' decreased by ' num2str(-deltaTs(i)*1e9) ' ns']);
            elseif deltaTs(i) > 0
                disp(['Duration of subinterval ' num2str(i) ' increased by ' num2str(deltaTs(i)*1e9) ' ns']);
            end
        end
    end
else
    disp('Some timing intervals were eliminated')
end

end