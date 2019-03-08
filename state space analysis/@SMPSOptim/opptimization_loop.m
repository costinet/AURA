function [] = opptimization_loop(obj)
%OPPTIMIZATION_LOOP will find the best time interval to achieve better
%losses and goal values of output parameters.
%   Detailed explanation goes here


%     _   _   _  ____    _
%    / \ | | | |/ _  |  / \
%   / _ \| | | | (_| | / _ \
%  / ___ | |_| |> _  |/ ___ \
% /_/   \_\___//_/ |_/_/   \_\

%% Constants from other classes

StateNumbers = obj.Simulator.Converter.Topology.Parser.StateNumbers;
StateNumbers_Opp = obj.Simulator.Converter.Topology.Parser.StateNumbers_Opposite;

Voerr = 999999999;
%% Adjust Power times

while Voerr>0.05*obj.Vo_ideal_value
    
    [y] = obj.Simulator.bruteforcelsim(100);
    
    ts = obj.Simulator.ts;
    Xs = obj.Simulator.Xs;
    As = obj.Simulator.As;
    Bs = obj.Simulator.Bs;
    u = obj.Simulator.u;
    
    
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
    
    
    
    % Cubic convergence
    %{
ts_temp = ts;
ts_new = ts;
dvodt = (mean(dXs(Vo_index,:))-mean(Xs(Vo_index,:)))/delta_DTs;
f = mean(Xs(Vo_index,:));
avg = f-Vo_ideal_value;
Newt = avg/dvodt;

ts_temp(Perturb1_index)=ts(Perturb1_index)-Newt;
ts_temp(Perturb2_index)=ts(Perturb2_index)+Newt;

[f_1int] = obj.Simulator.SS_Soln(0,As,Bs,ts_temp,u); % Zero is for keeping Xss from last time, not sure if this is applicable here. might want to pass a variable in the future
[f_1int] = obj.Simulator.CorrectXs(0,f_1int);

f2 = mean(f_1int(Vo_index,:));

avg2 = f2-Vo_ideal_value;

ts_new(Perturb1_index)=ts(Perturb1_index)-(Newt*avg)/(avg-avg2);
ts_new(Perturb2_index)=ts(Perturb2_index)+(Newt*avg)/(avg-avg2);
    %}
    
    %PHI = ts(1)-(u*f)/(f-)
    
    % Implement the change
    ts(1) = ts(1)+dt;
    ts(3) = ts(3)-dt;
    
    
    
end

end

