function describeInsertedIntervals(obj, allChanges)
% DESCRIBEINSERTEDINTERVALS is used as a debug function to describe if any
% intervals were inserted using findValidSteadyState when debug2=true
%
%   Example: 
%
%   If no changes were made the function will print to the command window:
%   No changes in switching pattern required, only times need to be adjusted
%
%   If there is a cHange made it will follow the format of: 
%   A new interval was inserted at the [end/beginning] of interval [interval
%   number] with [switch name] turned [on/off]
%
%   See also FINDVALIDSTEADYSTATE.




    if (isempty(allChanges))
        disp(['No changes in switching pattern required, only times need to be adjusted'])
    else
        interval = allChanges(:,1);
        beforeAfter = allChanges(:,2); 
        newSwInd = allChanges(:,3);
        switches = allChanges(:,4); 
        newStates = allChanges(:,5);
        for i = 1:size(allChanges,1) 
            if beforeAfter(i) == 1
                locString = 'at the end of';
            else
                locString = 'at the beginning of';
            end
    
            switchName = obj.switchNames{switches(i)};
    
            if newStates(i) == 1
                switchWord = 'on';
            else
                switchWord = 'off';
            end
    
            disp(['A new interval was inserted ' locString ' interval ' num2str(interval(i)) ' with ' switchName ' turned ' switchWord])
        end
    end
end