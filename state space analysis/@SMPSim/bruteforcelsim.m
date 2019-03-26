function [y] = bruteforcelsim(obj,iterations)
%BRUTEFORCELSIM Attempts to optimize the time intervals of a converter using lsim()
%   iterations is the number of iterations that will be looped through to try and optimize the converter provided in the class
%
%     _   _   _  ____    _    
%    / \ | | | |/ _  |  / \   
%   / _ \| | | | (_| | / _ \  
%  / ___ | |_| |> _  |/ ___ \ 
% /_/   \_\___//_/ |_/_/   \_\



% VDmax is the limited forward voltage allowed on a diode
% IDmax is the reverse current allowed on a diode
% VMmax is the limited forward voltage allowed on a diode
% IMmax is the reverse current allowed on a body diode (FET OFF)
% VDmax = 3;
% IDmax = -0.1;
% VMmax = -3;
% IMmax = 0.1;

tol = 0.005;

try
%tic
history_i = [];
history_j = [];
the_big_counter = 0;
not_reached_SS = true;
more_iterations=iterations;
while not_reached_SS && the_big_counter<=more_iterations
    the_big_counter = the_big_counter+1;
    goal_SS = obj.Xs(:,1);
    
    i = [];
    
    the_counter = 0;
    Xs = obj.Xs;
    ts = obj.ts;
    Ts = sum(ts);
    delta_t = 0;
    
    not_physical = true;
    
    while not_physical == 1 && the_counter<=iterations
        first_violations = false;
        last_violations = false;
        % pause(1)
        not_physical = false;
        the_counter = the_counter+1;
        Xs = obj.Xs;
        ts = obj.ts;
        Ts = sum(ts);
        % Important set up stuff
        debug = true;
        [xs, t, y, time_interval] = obj.SS_WF_Reconstruct();
        StateNumbers = obj.Converter.Topology.Parser.StateNumbers;
        StateNumbers_Opp = obj.Converter.Topology.Parser.StateNumbers_Opposite;
        ONorOFF = obj.Converter.Topology.Parser.ONorOFF;
        
        if debug
            fprintf('--------- \n')
            fprintf('Iteration number %.0f \n',the_counter)
            % What we are given for this iteration
            figure(1)
            ns = size(xs,1);
            for z=1:ns
                ax = subplot(10*ns,1,z*10-9:z*10);
                hold on;
                plot(t,y(StateNumbers(z),:), 'Linewidth', 3);
                ylabel(obj.getstatenames{z})
                box on
                ax.YLim = [min(y(StateNumbers(z),:))-abs(0.5*min(y(StateNumbers(z),:))) max(y(StateNumbers(z),:))+abs(0.5*max(y(StateNumbers(z),:)))];
                if(z<ns)
                    set(gca, 'Xticklabel', []);
                else
                    xlabel('t(s)')
                end
            end
            drawnow;
        end
        
        
        
       
        order = obj.order;
        time_ratio = Ts/time_interval(end); % Find the step value of the lsim that was done
        j = 2;
        new_index=1;
        time_variable_size = size(Xs,2);
        while j <= time_variable_size
            
            
            %{
        %% Universal Constraints
        if ONorOFF(i,1) ~=0 % if FET or Diode
            if obj.Converter.Topology.Parser.DMpos(i,2)==1 % if diode
                % Can never have voltage greater than forward voltage
                if ~isempty(find(y(StateNumbers(i),:) > VDmax, 1))&&debug % Only finds first violation to improve speed
                    fprintf('Universal Violation of diode forward voltage %s exceed %.2f V \n',obj.Converter.Topology.Parser.StateNames{i,1},VDmax)
                end
                % Can never have reverse current
                if ~isempty(find(y(StateNumbers_Opp(i),:) < IDmax,1)) && debug % Only finds first violation to improve speed
                    fprintf('Universal Violation of diode reverse current %s exceed %.2f A \n',obj.Converter.Topology.Parser.StateNames{i,1},IDmax)
                end
            end
        elseif obj.Converter.Topology.Parser.DMpos(i,3)==1 % if FET
            % Can never have voltage greater than forward voltage
            if ~isempty(find(y(StateNumbers(i),:) < VMmax, 1))&&debug % Only finds first violation to improve speed
                fprintf('Universal Violation of FET forward voltage %s exceed %.2f V \n',obj.Converter.Topology.Parser.StateNames{i,1},VMmax)
            end
            % Can have reverse current when FET is on will have to
            % check when FET is off individually by time interval
        end
            %}
            %% Cycle through time intervals
            

            for i = 1:1:size(Xs,1) % Cycle through state variables
                j;% is time interval for Xss
                k = j-1; % k is time interval for everything else
                
                if ONorOFF(i,k) ~=0 % if FET or Diode
                    if k==1
                        index = time_interval(k):time_interval(k+1); % List all the index values within the givn deadtime
                    else
                        index = time_interval(k)+1:time_interval(k+1); % List all the index values within the givn deadtime
                    end
                    waveform = y(StateNumbers(i),index); % Round off the values of the lsim
                    
                    if (time_ratio*length(index))~=ts(k)
                        somethingcool = ts(k)-(time_ratio*length(index));
                        delta_t = delta_t+ ts(k)-(time_ratio*length(index));
                    end
                    
                    if obj.Converter.Topology.Parser.DMpos(i,2)==1 % if diode
                        % Determine if a state changes from the Diode
                        % being on to being off or vice vera.
                        if ONorOFF(i,k) == 1 % if diode ON
                            if sum(waveform<1-1*tol)>0 && debug
                                fprintf('State Violation (Diode turn off) of %s in time interval %.0f \n',obj.Converter.Topology.Parser.StateNames{i,1},j-1)
                                if waveform(1)<1
                                    first_violations = 1;
                                    
                                    if the_counter>=5
                                        if  isequal(history_i(end),history_i(end-1),history_i(end-2),history_i(end-3),i) && isequal(history_j(end),history_j(end-1),history_j(end-2),history_j(end-3),j)
                                            %  old_ts_2 = obj.ts;
                                            if first_violations && ONorOFF(i,k-1) == -1
                                                j = j-1;
                                                [ts] = obj.Baxter_adjustDiodeConduction(obj.Xs,j,i,1,min(y(StateNumbers(i),:)));
                                                order = obj.order;
                                                not_physical = true;
                                                break
                                            end
                                        end
                                    end
                                    
                                end
                                sign = [0,0];
                                [ts,order,new_index] = obj.adjust_time(sign,new_index,ts,order,waveform,i,k,time_ratio);
                                not_physical = true;
                                break
                                %{
                            [P] = find(diff(waveform<1)~=0);
                            flipflop = waveform(1)<1;
                            ts_new = zeros(1,length(P)+1);
                            order_new = ts_new;
                            bd_state = obj.Converter.Topology.Parser.BD_OFF_state;
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
                            
                                %}
                                
                                %                             if all(diff(P)==1) && P(end)==size(waveform,2)
                                %                                 dt = time_ratio*(P(end)-P(1)+1);
                                %                                 obj.adjust_time(-ONorOFF(i,j-1),dt,k,i);
                                %                                 j = j+1;
                                %                             end
                            end
                            
                            
                        elseif ONorOFF(i,j-1) == -1 % if diode off
                            %  [V,P] = find(waveform > 1); % to try and find a value that is close to the goal deadtime value
                            if sum(waveform>1+1*tol)>0 && debug
                                fprintf('State Violation (Diode turn on) of %s in time interval %.0f \n',obj.Converter.Topology.Parser.StateNames{i,1},j-1)
                                
                                if waveform(end)>1
                                    last_violations = 1;
                                    
                                    
                                    if the_counter>=5
                                        if  isequal(history_i(end),history_i(end-1),history_i(end-2),history_i(end-3),i) && isequal(history_j(end),history_j(end-1),history_j(end-2),history_j(end-3),j)
                                            old_ts_2 = obj.ts;
                                            if last_violations && ONorOFF(i,k+1) == 1
                                                
                                                [ts] = obj.Baxter_adjustDiodeConduction(obj.Xs,j,i,1,min(y(StateNumbers(i),:)));
                                                order = obj.order;
                                                not_physical = true;
                                                break
                                            end
                                        end
                                    end
                                    
                                end
                                sign = [1,0];
                                [ts,order,new_index] = obj.adjust_time(sign,new_index,ts,order,waveform,i,k,time_ratio);
                                not_physical = true;
                                
                                break
                                %{
                            % Code to copy
                            [P] = find(diff(waveform>1)~=0);
                            flipflop = waveform(1)>1;
                            ts_new = zeros(1,length(P)+1);
                            order_new = ts_new;
                            bd_state = obj.Converter.Topology.Parser.BD_state;
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
                            
                                %}
                                
                                
                                %{
                            ts = obj.ts;
                            order = obj.order;
                            bd_state = obj.Converter.Topology.Parser.BD_state;
                            bd_off_state = obj.Converter.Topology.Parser.BD_OFF_state;
                            switches = obj.Converter.Topology.Parser.Switches;
                            ordered_state_index = obj.Converter.Topology.Parser.OrderedNamesnum;
                            
                            
                            bd_state_new = bd_state(order(time_index),ordered_state_index(k,order(time_index))==switches);
                            
                            if new_state == 1 % Turn body diode ON in a state
                                
                                
                                ts
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
                                
                                
                                %                             if all(diff(P)==1) && P(end)==size(waveform,2)
                                %                                 dt = time_ratio*(P(end)-P(1)+1);
                                %                                 obj.adjust_time(-ONorOFF(i,j-1),dt,k,i);
                                %                                 j = j+1;
                                %                             end
                            
                            end
                                %}
                            end
                        else
                            fprintf('Messed up\n')
                        end
                    elseif obj.Converter.Topology.Parser.DMpos(i,3)==1 % if FET
                        
                        
                        
                        if ONorOFF(i,j-1) == 2 % if FET ON
                            if sum(waveform<-1-1*tol)>0 && debug
                                % Check if body diode conducts (time interval could be shortened)
                                fprintf('Body Diode conducting: %s in time interval %.0f \n',obj.Converter.Topology.Parser.StateNames{i,1},j-1)
                            end
                            
                        elseif ONorOFF(i,j-1) == -1 % if FET off
                            
                            
                            
                            if sum(waveform<-1-1*tol)>0 && debug
                                
                                fprintf('State Violation (Body Diode turn on) of %s in time interval %.0f \n',obj.Converter.Topology.Parser.StateNames{i,1},j-1)
                                sign = [0,1];
                                [ts,order,new_index] = obj.adjust_time(sign,new_index,ts,order,waveform,i,k,time_ratio);
                                
                                not_physical = true;
                                if waveform(1)<-1
                                    first_violations = 1;
                                end
                                break
                                
                                %                             if all(diff(P)==1) && P(end)==size(waveform,2)
                                %                                 dt = time_ratio*(P(end)-P(1)+1);
                                %                                 obj.adjust_time(-ONorOFF(i,j-1),dt,k,i);
                                %                                 j = j+1;
                                %                             end
                                
                            end
                            
                            
                            
                        elseif ONorOFF(i,j-1) == 1 % body diode on
                            
                            if sum(waveform>-1+1*tol)>0 && debug
                                fprintf('Body Diode conducting: %s in time interval %.0f \n',obj.Converter.Topology.Parser.StateNames{i,1},j-1)
                                sign = [1,1];
                                [ts,order,new_index] = obj.adjust_time(sign,new_index,ts,order,waveform,i,k,time_ratio);
                                not_physical = true;
                                if waveform(1)>-1
                                    first_violations = 1;
                                end
                                break
                            end
                            
                            
                            %{
                        % Can not have positive ID
                        if ~isempty(find(y(StateNumbers_Opp(i),time_interval(i)+1:time_interval(i+1)) > IMmax,1)) && debug % Only finds first violation to improve speed
                            fprintf('State Violation of FET reverse current %s exceed %.2f A \n',obj.Converter.Topology.Parser.StateNames{i,1},IDmax)
                        end
                        
                        if ~isempty(find(y(StateNumbers_Opp(i),time_interval(i)+1:time_interval(i+1)) < -1,1)) && debug
                            fprintf('Body Diode conducting: %s in time interval %.0f \n',obj.Converter.Topology.Parser.StateNames{i,1},j-1)
                        end
                            %}
                            %                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %                         % Find if FET turns on next
                            %                         if size(ONorOFF,2)+1 == j
                            %                             if ONorOFF(i,2) == 1
                            %                                 if Voltage < -10*ron*2
                            %                                     fprintf('Hard swithcing for %s in time interval j\n',obj.Converter.Topology.Parser.StateNames{i,1},j-1)
                            %                                 end
                            %                             end
                            %
                            %                         else
                            %                             if ONorOFF(i,j) == 1
                            %                                 if Voltage < -10*ron*2
                            %                                     fprintf('Hard swithcing for %s in time interval j\n',obj.Converter.Topology.Parser.StateNames{i,1},j-1)
                            %                                 end
                            %                             end
                            %                         end
                            %                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        else
                            fprintf('Messed up\n')
                        end
                    else
                        fprintf('Found a FET or Diode that wasn''t a FET or Diode\n Don''t panic\n')
                    end
                    
                end

            end
            if not_physical
                break
            end
            j = j+1;
            new_index=new_index+1;
        end
        
        
        
        if not_physical % If there were no changes made then skip this step
            % Combine similar adjacent states would be nice!
            
            % This fixes the time intervals so there are not two adjacent time
            % intervals with the same state
            the_size = length(order);
            the_key = 1;
            while the_key < the_size
                if order(the_key)==order(the_key+1)
                    order(the_key+1) = [];
                    ts(the_key) = ts(the_key) + ts(the_key+1);
                    ts(the_key+1) = [];
                    the_size = the_size-1;
                    the_key = the_key-1;
                end
                the_key = the_key +1;
            end
            
            
            
            obj.ts = ts;
            obj.order = order;
            
            
            
            
            obj.updateTestConverter();
            obj.SS_Soln();
            obj.CorrectXs();
            obj.Converter.Topology.Parser.find_diode(obj.order);
        end
        
        history_i(end+1)=i;
        history_j(end+1)=j;
        
        %{
        if the_counter>=4
            if  isequal(history_i(end),history_i(end-1),history_i(end-2),i) && isequal(history_j(end),history_j(end-1),history_j(end-2),j)
               old_ts_2 = obj.ts;
               if first_violations
                   j = j-1;
               end
               [obj.ts] = obj.Baxter_adjustDiodeConduction(obj.Xs,j,i,max(y(StateNumbers(i),:)),min(y(StateNumbers(i),:)));
            end
        end
        history_i(end+1)=i;
        history_j(end+1)=j;
        
        obj.SS_Soln();
        obj.CorrectXs();
        obj.Converter.Topology.Parser.find_diode(obj.order);
        %}
    end
    check = 1;
    %obj.Xs(:,1) = obj.Xs(:,end); % Test step twords ss soln
    obj.SS_Soln(1);
    obj.CorrectXs(1);
    obj.Converter.Topology.Parser.find_diode(obj.order);
    
    not_reached_SS = ~isequal(goal_SS,obj.Xs(:,1));
    %obj.Baxter_adjustDiodeConduction(Xs,3,2,1,-100)
end

catch ME
    
    rethrow(ME)
end
%toc
end % That's all Folks



% Determine if it its possible to reach ideal value for soft switching (lsim for length of period)
% If it is possible to reach value then set the time to the first instance of that value
% If is not possible to reach value then step in direction by a certain percentage of initial guess.

%{
    
    if the_counter == 5
        J=234;
    end
    
    
    if mod(the_counter,2)
        %% Adjust Deadtimes
        deadtimes = obj.dead_time_intervals ;
        deadstates = obj.dead_time_states;
        deadgoals = obj.dead_time_goals;
        ts = obj.ts;
        Ts = sum(ts);
        
        % How do I translate between time in the SS_WF domain to the acutal
        % time lengths of each state that I need to set to run SSsoln.
        
        for i = 1:1:length(deadtimes)
            time_ratio = Ts/time_interval(end); % Find the step value of the lsim that was done
            index = time_interval(deadtimes(i))+1:time_interval(deadtimes(i)+1); % List all the index values within the givn deadtime
            waveform = round(y(StateNumbers(deadstates(i)),index),0); % Round off the values of the lsim
            [V,P] = find(waveform == deadgoals(i)); % to try and find a value that is close to the goal deadtime value
            
            % If it is not already fairly close
            if ~isempty(V) % If there was a match close to the deadtime goal

                if size(obj.As,3)==deadtimes(i)
                    ts(1) = ts(1)+ts(deadtimes(i))-time_ratio*P(1); % Set the next time interval to be longer accounting for the adjusted deadtime (special exception when deadtime is last time interval in period)
                else
                    ts(deadtimes(i)+1) = ts(deadtimes(i)+1)+ts(deadtimes(i))-time_ratio*P(1); % Set the next time interval to be longer accounting for the adjusted deadtime
                end
                ts(deadtimes(i))=time_ratio*P(1); % Set the deadtime to the found value (multiply index times the time ratio)
                
                % Set time for this interval to be this value if overshoot
            else

                % Expand seach to a period length
                As = obj.As; % get all the states from class
                Bs = obj.Bs;
                Cs = obj.Cs;
                Ds = obj.Ds;
                u = obj.u;
                test_t=sum(ts); % sum ts to find length of period
                ti = linspace(0,test_t,10000); % form the linspace for lsim solve
                SS = ss(As(:,:,deadstates(i)), Bs(:,:,deadstates(i)), Cs(:,:,deadstates(i)), Ds(:,:,deadstates(i))); % set up SS
                [y1, ~, ~] = lsim(SS, u*ones(size(ti)), ti, obj.Xs(:,deadtimes(i))); % run lsim
                time_ratio1 = Ts/10000; % find new Ratio of time to index values
                y1 = [y1'];
                % Same as steps above (except now with a longer time
                % interval), if the waveform is close to reaching the
                % goal then set it to it
                waveform = round(y1(StateNumbers(deadstates(i)),:),1);
                [V,P] = find(waveform == deadgoals(i));
                if ~isempty(V)
                    
                    if size(obj.As,3)==deadtimes(i)
                        ts(1) = ts(1)+ts(deadtimes(i))-time_ratio1*P(1);
                    else
                        ts(deadtimes(i)+1) = ts(deadtimes(i)+1)+ts(deadtimes(i))-time_ratio1*P(1);
                    end
                    ts(deadtimes(i))=time_ratio1*P(1);
                    % Set time for this interval to be this value if overshoot
                else
                    % If the value is not found then minimize error
                    % from goal
                    [V,P]=min(abs(deadgoals(i)-round(waveform,3)));

                    if size(obj.As,3)==deadtimes(i)
                        ts(1) = ts(1)+ts(deadtimes(i))-time_ratio1*P(1);
                    else
                        ts(deadtimes(i)+1) = ts(deadtimes(i)+1)+ts(deadtimes(i))-time_ratio1*P(1);
                    end
                    ts(deadtimes(i))=time_ratio*P(1);
                    % Move towards the region of less error (reduce error)
                end
                
            end
            
        end
        
        
    else
         %% Adjust Power times
        
        Vo_index = obj.Vo_index;
        Vo_ideal_value = obj.Vo_ideal_value;
        average = mean(y(StateNumbers(Vo_index),:)); % average output voltage
        Voerr = Vo_ideal_value-average; % Output voltage error
        Perturb1_index = obj.Perturb1_index;
        Perturb2_index = obj.Perturb2_index;
        
        % Preterb and observe approach using StateSensitivity.m
        delta_DTs = max(min(ts)/10, sum(ts)/10000);
        [dXs]=obj.Baxter_StateSensitivity('ts', Perturb1_index, delta_DTs,Perturb2_index); % Determine how much states change when time is changes
        dxsdt = (dXs-obj.Xs)/delta_DTs; % Find the change in states over the change in time (linear approximation of how changing ts will affect SS Xs)
        dt = (Voerr)/mean(dxsdt(Vo_index,:)); % Find the appropriate dt for output error
        
        % Implement the change
        ts(1) = ts(1)+(iterations-the_counter)/iterations*dt;
        ts(3) = ts(3)-(iterations-the_counter)/iterations*dt;
        
    end
%}
% if average>Vo_ideal_value
%     ts(1) = ts(1)-0.2*1/the_counter*ts(1);
%     ts(3) = ts(3)+0.2*1/the_counter*ts(1);
% end
%
%
% if average<Vo_ideal_value
%     ts(1) = ts(1)+0.2*1/the_counter*ts(1);
%     ts(3) = ts(3)-0.2*1/the_counter*ts(1);
% end

%obj.ts = ts;

