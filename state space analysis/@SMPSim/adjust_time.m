function [ts] = adjust_time(obj,new_state,dt,time_index,k,i,old_time)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

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
end

