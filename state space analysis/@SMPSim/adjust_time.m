function [ts,order,new_index] = adjust_time(obj,sign,new_index,ts,order,waveform,i,k,time_ratio)
%adjaust_time calcualtes the time needed and ordered states needed to
%adjust to have a correct physical circuit with diode implementation
%   Detailed explanation goes here



if sign(1)  == 1
    [P] = find(diff(waveform>1*(-1)^sign(2))~=0); % Find all cases of violations
    flipflop = waveform(1)>1*(-1)^sign(2); % Check to see if the inital value was in violation. This will be used later on if there are multiple crossings of violations within the same time interval. We simply flipflop between states in the order variable
    ts_new = zeros(1,length(P)+1); % intialize new ts two be inserted for the time 
    order_new = ts_new;
    % Get all of the values we need from Parser
    if sign(2) == 1
        bd_state = obj.Converter.Topology.Parser.BD_OFF_state;
    end
    if sign(2) == 0
        bd_state = obj.Converter.Topology.Parser.BD_state;
    end
    old_state = obj.order(k);
    ordered_state_index = obj.Converter.Topology.Parser.OrderedNamesnum;
    switches = obj.Converter.Topology.Parser.Switches;
    bd_state_new = bd_state(obj.order(k),ordered_state_index(i,obj.order(k))==switches);
    
    if isempty(P)
        % if there is no crossing then the entire state is in
        % violation of the
        order(new_index) = bd_state_new;
    else
        for number_of_changes = 1:1:length(P)+1
            order_new(number_of_changes) = flipflop;
            if number_of_changes == 1
                ts_new(1) = time_ratio*(P(1));
                
            elseif number_of_changes == length(P)+1
                ts_new(number_of_changes) = time_ratio*(length(waveform)-P(number_of_changes-1));
                
            else
                ts_new(number_of_changes) = time_ratio*(P(number_of_changes)-P(number_of_changes-1));
                
            end
            
            flipflop=~flipflop;
        end
        
        order_new(order_new==1) = bd_state_new(1);
        order_new(order_new==0) = old_state;
        ts = [ts(1:new_index-1) ts_new ts(new_index+1:end)];
        order = [order(1:new_index-1) order_new order(new_index+1:end)]; % Assume that FET is off and body diode is ON will have to adjust then when that correction is made its the (1) index

        new_index = new_index+length(P);
    end
end

if sign(1)  == 0
    [P] = find(diff(waveform<1*(-1)^sign(2))~=0);
    flipflop = waveform(1)<1*(-1)^sign(2);
    ts_new = zeros(1,length(P)+1);
    order_new = ts_new;
    
    if sign(2) == 0
        bd_state = obj.Converter.Topology.Parser.BD_OFF_state;
    end
    if sign(2) == 1
        bd_state = obj.Converter.Topology.Parser.BD_state;
    end
    
    old_state = obj.order(k);
    ordered_state_index = obj.Converter.Topology.Parser.OrderedNamesnum;
    switches = obj.Converter.Topology.Parser.Switches;
    bd_state_new = bd_state(obj.order(k),ordered_state_index(i,obj.order(k))==switches);
    if isempty(P)
        order(new_index) = bd_state_new;
    else
        for number_of_changes = 1:1:length(P)+1
            order_new(number_of_changes) = flipflop;
            if number_of_changes == 1
                ts_new(1) = time_ratio*(P(1));
                
            elseif number_of_changes == length(P)+1
                ts_new(number_of_changes) = time_ratio*(length(waveform)-P(number_of_changes-1));
                
            else
                ts_new(number_of_changes) = time_ratio*(P(number_of_changes)-P(number_of_changes-1));
                
            end
            
            flipflop=~flipflop;
        end
        
        order_new(order_new==1) = bd_state_new(1);
        order_new(order_new==0) = old_state;
        ts = [ts(1:new_index-1) ts_new ts(new_index+1:end)];
        order = [order(1:new_index-1) order_new order(new_index+1:end)]; % Assume that FET is off and body diode is ON will have to adjust then when that correction is made its the (1) index
        
        new_index = new_index+length(P);

    end
end









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
ts = obj.ts;
order = obj.order;
bd_state = obj.Converter.Topology.Parser.BD_state;
bd_off_state = obj.Converter.Topology.Parser.BD_OFF_state;
switches = obj.Converter.Topology.Parser.Switches;
ordered_state_index = obj.Converter.Topology.Parser.OrderedNamesnum;

if new_state == 1 % Turn body diode ON in a state
    
    ts(time_index) = ts(time_index)-dt;
    
    if ts(time_index)<5*eps(dt) % If you have essentially deleted the state
        ts(time_index) = 0;
        bd_state_new = bd_state(order(time_index),ordered_state_index(k,order(time_index))==switches);
        ts = [ts(1:time_index) dt ts(time_index+1:end)];
        order = [order(1:time_index) bd_state_new(1) order(time_index+1:end)];
        ts(time_index) = []; % Delete zero time states
        order(time_index) = []; % Delete zero time states
    else

    bd_state_new = bd_state(order(time_index),ordered_state_index(k,order(time_index))==switches);
    ts = [ts(1:time_index) dt ts(time_index+1:end)];
    order = [order(1:time_index) bd_state_new(1) order(time_index+1:end)]; % Assume that FET is off and body diode is ON will have to adjust then when that correction is made its the (1) index
    end
    
elseif new_state == -1 % Turn body diode OFF in a state
    
    ts(time_index) = ts(time_index)-dt;
    
    if ts(time_index)<5*eps(dt) % If you have essentially deleted the state
        ts(time_index) = 0;
        bd_off_state_new = bd_off_state(order(time_index),ordered_state_index(k,order(time_index))==switches);
        ts = [ts(1:time_index) dt ts(time_index+1:end)];
        order = [order(1:time_index) bd_off_state_new(1) order(time_index+1:end)];
        ts(time_index) = []; % Delete zero time states
        order(time_index) = []; % Delete zero time states
    else
        bd_off_state_new = bd_off_state(order(time_index),ordered_state_index(k,order(time_index))==switches);
        ts = [ts(1:time_index) dt ts(time_index+1:end)];
        order = [order(1:time_index) bd_off_state_new(1) order(time_index+1:end)]; % Assume that FET is off and body diode is ON will have to adjust then when that correction is made its the (1) index
    end
    
    
    

elseif newstate == 0 % Change the length of already existing states
    
    
    fprintf('Not done yet\n Be patient\n')
    
    
end

obj.order = order;
obj.ts = ts;

%}
end

