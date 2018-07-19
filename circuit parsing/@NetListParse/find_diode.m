function [] = find_diode(obj,order)
% Function to determine where diodes are in matrix
% this works!!

% DMpos is the output
% First column if 1 if it is a D or M
% Second column is 1 if D
% Thirt column is 1 if M

StateNames=obj.StateNames;
num=size(StateNames,1);

DMpos = zeros(num,3);

for i = 1:1:num % Iterate through number of states variables
    test=StateNames{i,1};
    if strcmp(test(1),'D') % if there is a match between the dependent name and the measurement name
        DMpos(i,:) = [1 1 0];
    elseif strcmp(test(1),'M')
        DMpos(i,:) = [1 0 1];
    else
        DMpos(i,:) = [0 0 0]; % Not neccessary
        % Component is not DM
    end
end

onoroff = zeros(num,length(order));

for i = 1:1:length(order) % Iterate through number of states
    for j = 1:1:num % Iterate through number of states variables
        for k = 1:1:size(obj.ON_States,1)
            if strcmp(obj.ON_States{k,order(i)},obj.StateNames(j,1))
                onoroff(j,i)=1;
            end
            if strcmp(obj.OFF_States{k,order(i)},obj.StateNames(j,1))
                onoroff(j,i)=-1;
            end
        end
    end
end


        
        


if sum(sum(DMpos,2)==1)~=0
    fprintf('Either someone messed with this code or all logic in the world has been lost\n Let''s hope its the first one\n')
end

obj.DMpos = DMpos;
obj.ONorOFF = onoroff;
end % That's all Folks

