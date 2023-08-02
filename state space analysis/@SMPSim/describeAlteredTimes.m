function describeAlteredTimes(obj,oldts)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

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