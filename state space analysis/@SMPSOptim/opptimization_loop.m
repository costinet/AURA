function [] = opptimization_loop(obj)
%OPPTIMIZATION_LOOP will find the best time interval to achieve better
%losses and goal values of output parameters.
%   Detailed explanation goes here

%% Constants from other classes

StateNumbers = obj.Simulator.Converter.Topology.Parser.StateNumbers;
StateNumbers_Opp = obj.Simulator.Converter.Topology.Parser.StateNumbers_Opposite;



%% Adjust Power times

[y] = obj.Simulator.bruteforcelsim(5);

ts = obj.Simulator.ts;
Xs = obj.Simulator.Xs;

Vo_index = obj.Vo_index;
Vo_ideal_value = obj.Vo_ideal_value;
average = mean(y(StateNumbers(Vo_index),:)); % average output voltage
Voerr = Vo_ideal_value-average; % Output voltage error
Perturb1_index = obj.Perturb1_index;
Perturb2_index = obj.Perturb2_index;

% Preterb and observe approach using StateSensitivity.m
delta_DTs = max(min(ts)/10, sum(ts)/10000);
[dXs]=obj.Simulator.Baxter_StateSensitivity('ts', Perturb1_index, delta_DTs,Perturb2_index); % Determine how much states change when time is changes
dxsdt = (dXs-Xs)/delta_DTs; % Find the change in states over the change in time (linear approximation of how changing ts will affect SS Xs)
dt = (Voerr)/mean(dxsdt(Vo_index,:)); % Find the appropriate dt for output error

% Implement the change
ts(1) = ts(1)+(iterations-the_counter)/iterations*dt;
ts(3) = ts(3)-(iterations-the_counter)/iterations*dt;




end

