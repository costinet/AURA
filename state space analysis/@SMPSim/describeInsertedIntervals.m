function describeInsertedIntervals(obj, allChanges)
%describeInsertedIntervals Print to console a description of how switching
%state order has changed
%   
%   describeInsertedIntervals(obj, allChanges)
%   prints out a description of how an SMPSim objects current switching
%   interval sequence was changed through the variable allChanges.  See
%   debug blocks in SMPSim.findValidSteadyState for correct usage
%
%   See Also SMPSim.describeDiscreteErrors,
%   SMPSim.describeAlteredTimes, SMPSim.describeSwitchState

    

    if (isempty(allChanges))
        disp(['No changes in switching pattern required, only times need to be adjusted'])
    else
        interval = allChanges(:,1);
        beforeAfter = allChanges(:,2); 
        newSwInd = allChanges(:,3);
        switches = allChanges(:,4); 
        newStates = allChanges(:,5);
        for loc = [unique([interval beforeAfter], 'rows')]'
            int = loc(1);
            BA = loc(2);

            rows = find(interval == int & beforeAfter == BA);
        
            if BA == 1
                locString = 'at the end of';
            else
                locString = 'at the beginning of';
            end

            switchString = '';
            for j = 1:length(rows)
                switchName = obj.switchNames{switches(rows(j))};
                if newStates(rows(j)) == 1
                    switchWord = 'on';
                else
                    switchWord = 'off';
                end
                    
                switchString = [switchString switchName ' turned ' switchWord];

                if j < length(rows)
                    switchString = [switchString ' and '];
                end
            end
    

    
            disp(['A new interval was inserted ' locString ' interval ' num2str(int) ' with ' switchString])
        end
    end
end