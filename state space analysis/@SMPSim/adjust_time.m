function [outputArg1,outputArg2] = adjust_time(obj,new_state,dt,time_index,k,i,old_time)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

ts = obj.ts;


if new_state == 1 % Turn body diode ON in a state

    
    dt == ts(k)


elseif newstate == -1 % Turn body diode OFF in a state
    
    
    
    
    

elseif newstate == 0 % Change the length of already existing states
    
    
    fprintf('Not done yet\n Be patient\n')
    
    
end

obj.ts
end

