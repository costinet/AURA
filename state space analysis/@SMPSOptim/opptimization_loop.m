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

%% Adjust Power times

    obj.Simulator.As_OG = obj.Simulator.As;
    obj.Simulator.Bs_OG = obj.Simulator.Bs;
    obj.Simulator.Cs_OG = obj.Simulator.Cs;
    obj.Simulator.Ds_OG = obj.Simulator.Ds;
    obj.Simulator.u_OG = obj.Simulator.u;
    obj.Simulator.eigA_OG = obj.Simulator.eigA;
    obj.Simulator.ONorOFF_OG = obj.Simulator.Converter.Topology.Parser.ONorOFF;
    obj.Simulator.ts_OG = obj.Simulator.ts;
    
onward = true;

while onward
    onward = false;
    obj.Simulator.Three_tier_diode_correct(40,0);
    
    ts = obj.Simulator.ts;
    Xs = obj.Simulator.Xs;
    As = obj.Simulator.As;
    Bs = obj.Simulator.Bs;
    u = obj.Simulator.u;
    ONorOFF = obj.Simulator.Converter.Topology.Parser.ONorOFF;
    
    
%     cumulativesum=cumsum(ts)./sum(ts);
%     findcumulative = find(abs(cumulativesum-0.5)<eps*100);
%     
    half_way_ONofOFF  =  [2 ;    2 ;   -1;     2  ;  -1  ;   0  ;  -1   ;  0   ;  0];
    
    
    ONorOFFtofind =  [ -1 ;    2  ;   2  ;  -1  ;  -1   ;  0  ;   2  ;   0  ;   0];
    
    indextochange= [];
    for i = 1:size(ONorOFF,2)
        if ONorOFF(:,i)== ONorOFFtofind
            indextochange(end+1) = i;
        end
        
        if ONorOFF(:,i)== half_way_ONofOFF
            findcumulative = i-1;
        end
            
    end
    
    Xsindex = 6; %This is to find the inductor row in the steady state solution
    GoaliL = 4.9; % This is the goal minimum inductor current
    Sir=Xsindex;
    

    
    for z = 1:2
        
        if z == 1
            Xic = findcumulative+1;
        elseif z==2
            Xic  = size(Xs,2);
        end
        
        
        
        if (abs(Xs(Xsindex,Xic)-GoaliL)>0.1) 
            onward = true;
            progBar = 0.001;
            maxStep = sum(ts)/100/(.5+1.5*progBar);
            
            
            
            
            Ti = indextochange(z);
            

            
            
            delta_DTs = max(min(ts)/1000, sum(ts)/100000);
            [dXs,delta_DTs] = obj.Simulator.Baxter_StateSensitivity2(0, 'ts', Ti, delta_DTs);
            dxsdt = (dXs-Xs)/delta_DTs;
            
            tdelta = (GoaliL-Xs(Sir,Xic))/(dxsdt(Sir,Xic));
            % tdelta = min(max(tdelta, -ts(Ti) + delta_DTs), tsmax(Ti) - ts(Ti));
            tdelta = sign(tdelta)*min(abs(tdelta), maxStep);
            ts(Ti) = ts(Ti) + tdelta;
            if(ts(Ti)<0)
                change = abs(ts(Ti))*0.1;
                ts(Ti) = change;
            end
            
            if z ==1
                obj.Simulator.ts_OG(5) = obj.Simulator.ts_OG(5)+tdelta;
            end
            
            
            if z ==2
                obj.Simulator.ts_OG(11) = obj.Simulator.ts_OG(11)+tdelta;
            end
            break
        end
    end
    
    obj.Simulator.setts(ts);
    obj.Simulator.SS_Soln();
    obj.Simulator.CorrectXs();
    
    
    
    
    
    
    
    
    
    %%%%% This is for the Buck Boost Converter only to do
    %%%%% switch at the current dependent point
%     if i==6 && abs((sum(ts(1:k)))-(sum(ts)/2))<eps*100
%         Hello = 5465;
%     end
    
    % Set switch inductor value to 5
    
    
    
    
    
    
    
    
    %[dXs] = obj.fs_adjust();
    
    
    % How to find the partial with respect to a variable 
    
    % For a range of Iout at 3 6 9 and 12 watts what is the partial with repect to Ron and
    % Roff at different 
    
    
    
    %{
    Vo_index = obj.Vo_index;
    Vo_ideal_value = obj.Vo_ideal_value;
    average = mean(y(StateNumbers(Vo_index),:)); % average output voltage
    Voerr = Vo_ideal_value-average; % Output voltage error
    Perturb1_index = obj.Perturb1_index;
    Perturb2_index = obj.Perturb2_index;
    
    % Preterb and observe approach using StateSensitivity.m
    delta_DTs = max(min(ts)/10, sum(ts)/10000);
    [dXs]=obj.Simulator.Baxter_StateSensitivity(0,'ts', Perturb1_index, delta_DTs,Perturb2_index); % Determine how much states change when time is changes
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
    ts(Perturb1_index) = ts(Perturb1_index)+dt;
    ts(Perturb2_index) = ts(Perturb2_index)-dt;
    
    obj.ts = ts;
   %} 
end
obj.Simulator.Three_tier_diode_correct(40,0);
end

