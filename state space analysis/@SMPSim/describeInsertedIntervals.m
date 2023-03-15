function describeInsertedIntervals(obj, allChanges)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    

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