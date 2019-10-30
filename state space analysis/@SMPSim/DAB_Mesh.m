function [] = DAB_Mesh(obj)
% DAB_Mesh takes the a DAB converter and assuming that the primary and
% switchign states are both symetric and want to be equivalent, then
% there can be a surface! The beginning of this code is simply used to
% place the deadtimes in the circuit fromt the parser, it is reused
% from bruteforecelsim_single

%
%     _   _   _  ____    _
%    / \ | | | |/ _  |  / \
%   / _ \| | | | (_| | / _ \
%  / ___ | |_| |> _  |/ ___ \
% /_/   \_\___//_/ |_/_/   \_\

iterations =1;
tol = 0.005;
toc
try
    %tic
    history_i = [];
    history_j = [];
    the_big_counter = 0;
    not_reached_SS = true;
    more_iterations=iterations;
    ONorOFF = obj.Converter.Topology.Parser.ONorOFF;
    
    the_big_counter = the_big_counter+1;
    goal_SS = obj.Xs(:,1);
    
    i = [];
    
    savedts =[];
    the_counter = 0;
    Xs = obj.Xs;
    ts = obj.ts;
    Ts = sum(ts);
    delta_t = 0;
    
    keep_SS = false;
    
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
        if debug
            output_order =  [10 8 6 9];
            [xs, t, y, time_interval] = obj.SS_WF_Reconstruct();
            StateNumbers = obj.Converter.Topology.Parser.StateNumbers;
            StateNumbers_Opp = obj.Converter.Topology.Parser.StateNumbers_Opposite;
        end
        ONorOFF = obj.Converter.Topology.Parser.ONorOFF;
        
        
        
        %%%
        % New variables to try and speed things up a little
        
        last_violations_bd_turn_on = zeros(size(ONorOFF));
        last_violations_bd_turn_off = last_violations_bd_turn_on;
        first_violations_bd_turn_off = last_violations_bd_turn_off;
        first_violations_bd_turn_on = last_violations_bd_turn_off;
        %%%
        
        %%% This is original waveform generator for all state
        %%% variables
        
        if debug
            fprintf('--------- \n')
            fprintf('Iteration number %.0f \n',the_counter-1)
            savedts(end+1,:) = ts;
        end
        if debug
            
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
        
        if debug
            % What we are given for this iteration
            figure(2)
            ns = size(xs,1);
            for z=1:ns
                ax = subplot(10*ns,1,z*10-9:z*10);
                hold on;
                plot(t,y(StateNumbers_Opp(z),:), 'Linewidth', 3);
                ylabel(obj.getstatenames_Opp{z})
                box on
                ax.YLim = [min(y(StateNumbers_Opp(z),:))-abs(0.5*min(y(StateNumbers_Opp(z),:))) max(y(StateNumbers_Opp(z),:))+abs(0.5*max(y(StateNumbers_Opp(z),:)))];
                if(z<ns)
                    set(gca, 'Xticklabel', []);
                else
                    xlabel('t(s)')
                end
            end
            drawnow;
        end
        
        
        %                 %%% This is new waveform generator for the listed
        %                 %%% waveforms of the compel converter
        %             if debug
        %                 fprintf('--------- \n')
        %                 fprintf('Iteration number %.0f \n',the_counter-1)
        %                 % What we are given for this iteration
        %                 figure(5)
        %                 ns = size(xs,1);
        %                 for z=1:length(output_order)
        %                     ax = subplot(10*length(output_order),1,z*10-9:z*10);
        %                     hold on;
        %                     plot(t,y(StateNumbers(output_order(z)),:), 'Linewidth', 3);
        %                     ylabel(obj.getstatenames{output_order(z)})
        %                     box on
        %                     ax.YLim = [min(y(StateNumbers(output_order(z)),:))-abs(0.5*min(y(StateNumbers(output_order(z)),:))) max(y(StateNumbers(output_order(z)),:))+abs(0.5*max(y(StateNumbers(output_order(z)),:)))];
        %                     if(z<length(output_order))
        %                         set(gca, 'Xticklabel', []);
        %                     else
        %                         xlabel('t(s)')
        %                     end
        %                 end
        %                 drawnow;
        %             end
        
        
        
        not_physical = true;
        the_little_counter = 1;
        
        %             while not_physical || the_little_counter>5
        the_little_counter = the_little_counter+1;
        not_physical = false;
        order = obj.order;
        if debug
            time_ratio = Ts/time_interval(end); % Find the step value of the lsim that was done
        end
        j = 2;
        time_variable_size = size(Xs,2);
        while j <= time_variable_size
            
            
            
            
            %% Cycle through time intervals
            
            
            for i = 1:1:size(Xs,1) % Cycle through state variables
                j;% is time interval for Xss
                k = j-1; % k is time interval for everything else
                
                if i ==9 && k==5
                    J = 45465;
                end
                
                
                if ONorOFF(i,k) ~=0 % if FET or Diode
                    if debug
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
                    end
                    if obj.Converter.Topology.Parser.DMpos(i,2)==1 % if diode
                        % Determine if a state changes from the Diode
                        % being on to being off or vice vera.
                        
                        if ONorOFF(i,k) == 1 % if diode ON
                            if sum(waveform<1-1*tol)>0  % if diode should turn off during time interval
                                if debug
                                    fprintf('State Violation (Diode turn off) of %s in time interval %.0f \n',obj.Converter.Topology.Parser.StateNames{i,1},j-1)
                                end
                                sign = [0,0];
                                [P] = find(diff(waveform<1*(-1)^sign(2))~=0); % Find all cases of violations
                                if length(P)==1||length(P)==2
                                    if waveform(1)<1 % if there is only a violation at the beginning of the time interval then look to set the end of the last time interval to be equal to the diode forward votlage (1)
                                        first_violations = 1;
                                        
                                        %                                     if the_counter>=5
                                        %                                         if  isequal(history_i(end),history_i(end-1),history_i(end-2),history_i(end-3),i) && isequal(history_j(end),history_j(end-1),history_j(end-2),history_j(end-3),j)
                                        %  old_ts_2 = obj.ts;
                                        if first_violations && ONorOFF(i,k-1) == -1
                                            j = j-1;
                                            [ ts, ~, ~, ~, ~,keep_SS] = obj.Baxter_adjustDiodeConduction(obj.Xs,j,i,1,min(y(StateNumbers(i),:)),0.001,1,keep_SS);
                                            order = obj.order;
                                            not_physical = true;
                                            break
                                        end
                                        %                                         end
                                        %                                     end
                                    end
                                    
                                    if waveform(end)<1 % if there is only a violation at the end of the time interval then the current time interval needs to be ajusted to end earilier at the diode forward votlage crossing (1)
                                        last_violations = 1;
                                        
                                        %                                     if the_counter>=5
                                        %                                         if  isequal(history_i(end),history_i(end-1),history_i(end-2),history_i(end-3),i) && isequal(history_j(end),history_j(end-1),history_j(end-2),history_j(end-3),j)
                                        %  old_ts_2 = obj.ts;
                                        if last_violations && ONorOFF(i,k+1) == -1
                                            [ ts, ~, ~, ~, ~,keep_SS] = obj.Baxter_adjustDiodeConduction(obj.Xs,j,i,max(y(StateNumbers(i),:)),1,0.001,0,keep_SS);
                                            order = obj.order;
                                            not_physical = true;
                                            break
                                        end
                                        %                                         end
                                        %                                     end
                                    end
                                    
                                else
                                    
                                    [ts,order,new_index] = obj.adjust_time(sign,new_index,ts,order,waveform,i,k,time_ratio);
                                    not_physical = true;
                                    break
                                end
                                
                            end
                            
                            
                        elseif ONorOFF(i,j-1) == -1 % if diode off
                            %  [V,P] = find(waveform > 1); % to try and find a value that is close to the goal deadtime value
                            if sum(waveform>1+1*tol)>0
                                if debug
                                    fprintf('State Violation (Diode turn on) of %s in time interval %.0f \n',obj.Converter.Topology.Parser.StateNames{i,1},j-1)
                                end
                                sign = [1,0];
                                [P] = find(diff(waveform>1*(-1)^sign(2))~=0); % Find all cases of violations
                                if length(P)==1
                                    if waveform(end)>1 % if the violation is at the end the need to reduce time interval to diode forward votlage crossing (1)
                                        last_violations = 1;
                                        
                                        
                                        %                                     if the_counter>=5
                                        %                                         if  isequal(history_i(end),history_i(end-1),history_i(end-2),history_i(end-3),i) && isequal(history_j(end),history_j(end-1),history_j(end-2),history_j(end-3),j)
                                        old_ts_2 = obj.ts;
                                        if last_violations && ONorOFF(i,k+1) == 1
                                            
                                            [ ts, ~, ~, ~, ~,keep_SS] = obj.Baxter_adjustDiodeConduction(obj.Xs,j,i,1,min(y(StateNumbers(i),:)),0.001,1,keep_SS);
                                            order = obj.order;
                                            not_physical = true;
                                            break
                                        end
                                        %                                         end
                                        %                                     end
                                        
                                    end
                                    
                                    [ts,order,new_index] = obj.adjust_time(sign,new_index,ts,order,waveform,i,k,time_ratio);
                                    not_physical = true;
                                    
                                    break
                                end
                            end
                        else
                            fprintf('Messed up\n')
                        end
                    elseif obj.Converter.Topology.Parser.DMpos(i,3)==1 % if FET
                        
                        
                        
                        if ONorOFF(i,j-1) == 2 % if FET ON
                            %                                 if sum(waveform<-1-1*tol)>0
                            %                                     Hello_world  = 9;
                            %                                     % Check if body diode conducts (time interval could be shortened)
                            %                                     %  fprintf('Body Diode conducting: %s in time interval %.0f \n',obj.Converter.Topology.Parser.StateNames{i,1},j-1)
                            %                                 end
                            
                        elseif ONorOFF(i,j-1) == -1 % if FET off
                            
                            
                            if (obj.Xs(i,j)<-1-1*tol || obj.Xs(i,j-1)<-1-1*tol )
                                
                                if debug
                                    fprintf('State Violation (Body Diode turn on) of %s in time interval %.0f \n',obj.Converter.Topology.Parser.StateNames{i,1},j-1)
                                end
                                sign = [0,1];
                                
                                
                                %                                     [P] = find(diff(waveform<1*(-1)^sign(2))~=0); % Find all cases of violations
                                
                                %                                     if (isempty(P))
                                %                                         % Then there has been a
                                %                                         % violation detected but no
                                %                                         % crossover points so the
                                %                                         % entire time interval is in
                                %                                         % error
                                %                                         replace = 1;
                                %                                         [ts,order,new_index] = obj.adjust_time_single(sign,new_index,ts,order,waveform,i,k,time_ratio,ONorOFF,replace);
                                %
                                %                                         not_physical = true;
                                %                                         fprintf('All of state was in error \n');
                                %
                                %                                         break
                                %                                     end
                                
                                
                                
                                %                                     if (length(P)<2)
                                
                                if obj.Xs(i,j)<-1 % if there is only a violation at the end of the time interval then the current time interval needs to be ajusted to end earilier at the diode forward votlage crossing (1)
                                    last_violations = 1;
                                    if (k+1)>size(ONorOFF,2)
                                        %         if last_violations %&& ONorOFF(i,1) == 2
                                        %           [ ts, ~, ~, ~, ~,keep_SS] = obj.Baxter_adjustDiodeConduction(obj.Xs,j,i,max(y(StateNumbers(i),:)),-1,0.001,0,keep_SS);
                                        %         order = obj.order;
                                        %       not_physical = true;
                                        %     break
                                        %      end
                                        
                                        if ONorOFF(i,1) == 1 && (sum((ONorOFF(:,1)==2)==(ONorOFF(:,j-1)==2)) == size(ONorOFF,1)) % This checks to see if there is already a diode state for the next interval that does not affect switching actions
                                            [ ts, ~, ~, ~, ~,keep_SS] = obj.Baxter_adjustDiodeConduction(obj.Xs,j,i,max(y(StateNumbers(i),:)),-1,0.001,0,keep_SS);
                                            obj.setts(ts);
                                            order = obj.order;
                                            not_physical = true;
                                            %   obj.SS_Soln(1);
                                            %  obj.CorrectXs(1);
                                            
                                            
                                        else
                                            %                                                     replace = 0; new_index = 1;
                                            %                                                     [ts,order,new_index] = obj.adjust_time_single(sign,new_index,ts,order,waveform,i,k,time_ratio,ONorOFF,replace);
                                            %                                                     ts = [ts(1:k-1) ts(k)*.95 ts(k)*.05 ts(k+1:end)];
                                            %                                                     obj.setts(ts);
                                            %                                                     obj.SS_Soln(keep_SS);
                                            %                                                     obj.CorrectXs(keep_SS);
                                            %                                                     [ ts, ~, ~, ~, ~,keep_SS] = obj.Baxter_adjustDiodeConduction(obj.Xs,j,i,max(y(StateNumbers(i),:)),-1,0.001,0,keep_SS);
                                            %                                                     obj.setts(ts);
                                            %                                                     order = obj.order;
                                            %                                                     not_physical = true;
                                            %                                                     break
                                            
                                            last_violations_bd_turn_on(i,j-1) = 1;
                                            
                                        end
                                        
                                        
                                    else
                                        if last_violations %&& ONorOFF(i,k+1) == 2
                                            if ONorOFF(i,j) == 1 && (sum((ONorOFF(:,j)==2)==(ONorOFF(:,j-1)==2)) == size(ONorOFF,1)) % This checks to see if there is already a diode state for the next interval that does not affect switching actions
                                                [ ts, ~, ~, ~, ~,keep_SS] = obj.Baxter_adjustDiodeConduction(obj.Xs,j,i,max(obj.Xs(i,:)),-1,0.001,0,keep_SS);
                                                obj.setts(ts);
                                                order = obj.order;
                                                not_physical = true;
                                                %    obj.SS_Soln(1);
                                                %  obj.CorrectXs(1);
                                                
                                            else
                                                %                                                         replace = 0; new_index = 1;
                                                %                                                         [ts,order,new_index] = obj.adjust_time_single(sign,new_index,ts,order,waveform,i,k,time_ratio,ONorOFF,replace);
                                                %                                                         ts = [ts(1:k-1) ts(k)*.95 ts(k)*.05 ts(k+1:end)];
                                                %                                                         obj.setts(ts);
                                                %                                                         obj.SS_Soln(keep_SS);
                                                %                                                         obj.CorrectXs(keep_SS);
                                                %                                                         [ ts, ~, ~, ~, ~,keep_SS] = obj.Baxter_adjustDiodeConduction(obj.Xs,j,i,max(y(StateNumbers(i),:)),-1,0.001,0,keep_SS);
                                                %                                                         obj.setts(ts);
                                                %                                                         order = obj.order;
                                                %                                                         not_physical = true;
                                                %                                                         break
                                                
                                                last_violations_bd_turn_on(i,j-1) = 1;
                                                
                                            end
                                        end
                                    end
                                end
                                
                                if obj.Xs(i,j-1)<-1 % if there is only a violation at the beginning of the time interval then look to set the end of the last time interval to be equal to the diode forward votlage (1)
                                    first_violations = 1;
                                    
                                    if j == 2
                                        
                                        if first_violations %&& %ONorOFF(i,end) == 2
                                            
                                            if ONorOFF(i,end) == 1 && (sum((ONorOFF(:,j-1)==2)==(ONorOFF(:,end)==2)) == size(ONorOFF,1)) % This checks to see if there is already a diode state for the next interval that does not affect switching actions
                                                [ ts, ~, ~, ~, ~,keep_SS] = obj.Baxter_adjustDiodeConduction(obj.Xs,j-1,i,-1,min(obj.Xs(i,:)),0.001,1,keep_SS);
                                                obj.setts(ts);
                                                order = obj.order;
                                                not_physical = true;
                                                % obj.SS_Soln(1);
                                                %  obj.CorrectXs(1);
                                                
                                            else
                                                %                                                         replace = 0; new_index = 0;
                                                %                                                         [ts,order,new_index] = obj.adjust_time_single(sign,new_index,ts,order,waveform,i,k,time_ratio,ONorOFF,replace);
                                                %                                                         ts = [ts(1:k-1) ts(k)*.05 ts(k)*.95 ts(k+1:end)];
                                                %                                                         obj.setts(ts);
                                                %                                                         obj.SS_Soln(keep_SS);
                                                %                                                         obj.CorrectXs(keep_SS);
                                                %                                                         [ ts, ~, ~, ~, ~,keep_SS] = obj.Baxter_adjustDiodeConduction(obj.Xs,j,i,-1,min(y(StateNumbers(i),:)),0.001,1,keep_SS);
                                                %                                                         obj.setts(ts);
                                                %                                                         order = obj.order;
                                                %                                                         not_physical = true;
                                                %                                                         break
                                                
                                                first_violations_bd_turn_on(i,j-1) = 1;
                                                
                                            end
                                            
                                        end
                                        
                                    else
                                        
                                        if first_violations %&& %ONorOFF(i,k-1) == 2
                                            
                                            
                                            if ONorOFF(i,j-2) == 1 && (sum((ONorOFF(:,j-1)==2)==(ONorOFF(:,j-2)==2)) == size(ONorOFF,1)) % This checks to see if there is already a diode state for the next interval that does not affect switching actions
                                                [ ts, ~, ~, ~, ~,keep_SS] = obj.Baxter_adjustDiodeConduction(obj.Xs,j-1,i,-1,min(obj.Xs(i,:)),0.001,1,keep_SS);
                                                obj.setts(ts);
                                                order = obj.order;
                                                not_physical = true;
                                                % obj.SS_Soln(1);
                                                %  obj.CorrectXs(1);
                                                
                                            else
                                                %                                                         replace = 0; new_index = 0;
                                                %                                                         [ts,order,new_index] = obj.adjust_time_single(sign,new_index,ts,order,waveform,i,k,time_ratio,ONorOFF,replace);
                                                %                                                         ts = [ts(1:k-1) ts(k)*.05 ts(k)*.95 ts(k+1:end)];
                                                %                                                         obj.setts(ts);
                                                %                                                         obj.SS_Soln(keep_SS);
                                                %                                                         obj.CorrectXs(keep_SS);
                                                %                                                         [ ts, ~, ~, ~, ~,keep_SS] = obj.Baxter_adjustDiodeConduction(obj.Xs,j,i,-1,min(y(StateNumbers(i),:)),0.001,1,keep_SS);
                                                %                                                         obj.setts(ts);
                                                %                                                         order = obj.order;
                                                %                                                         not_physical = true;
                                                %                                                         break
                                                
                                                first_violations_bd_turn_on(i,j-1) = 1;
                                                
                                            end
                                            
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                                
                                
                                %                                     else
                                %
                                %
                                %
                                %
                                %
                                %                                         [ts,order,new_index] = obj.adjust_time(sign,new_index,ts,order,waveform,i,k,time_ratio);
                                %
                                %                                         not_physical = true;
                                %                                     end
                                
                                
                                %                             if all(diff(P)==1) && P(end)==size(waveform,2)
                                %                                 dt = time_ratio*(P(end)-P(1)+1);
                                %                                 obj.adjust_time(-ONorOFF(i,j-1),dt,k,i);
                                %                                 j = j+1;
                                %                             end
                                
                            end
                            
                            
                            
                        elseif ONorOFF(i,j-1) == 1 % body diode on
                            
                            if (obj.Xs(i,j)>-1+1*tol || obj.Xs(i,j-1)>-1+1*tol)
                                if debug
                                    fprintf('State Violation (Body Diode turn off) of %s in time interval %.0f \n',obj.Converter.Topology.Parser.StateNames{i,1},j-1)
                                end
                                
                                sign = [1,1];
                                
                                
                                %  [P] = find(diff(waveform>1*(-1)^sign(2))~=0); % Find all cases of violations
                                
                                %                                     if (isempty(P))
                                %                                         % Then there has been a
                                %                                         % violation detected but no
                                %                                         % crossover points so the
                                %                                         % entire time interval is in
                                %                                         % error
                                %                                         replace = 1;
                                %                                         [ts,order,new_index] = obj.adjust_time_single(sign,new_index,ts,order,waveform,i,k,time_ratio,ONorOFF,replace);
                                %
                                %                                         not_physical = true;
                                %                                         fprintf('All of state was in error \n');
                                %
                                %                                         break
                                %                                     end
                                %
                                
                                
                                %                                     if (length(P)<2)
                                
                                if obj.Xs(i,j)>-1 % if there is only a violation at the end of the time interval then the current time interval needs to be ajusted to end earilier at the diode forward votlage crossing (1)
                                    last_violations = 1;
                                    if (k+1)>size(ONorOFF,2)
                                        if last_violations %&& ONorOFF(i,1) == 2
                                            [ ts, ~, ~, ~, ~,keep_SS] = obj.Baxter_adjustDiodeConduction(obj.Xs,j,i,max(obj.Xs(i,:)),-1,0.001,0,keep_SS);
                                            order = obj.order;
                                            not_physical = true;
                                            % obj.SS_Soln(1);
                                            %       obj.CorrectXs(1);
                                        end
                                    else
                                        if last_violations %&& ONorOFF(i,k+1) == 2
                                            if ONorOFF(i,j) == -1 && (sum((ONorOFF(:,j)==2)==(ONorOFF(:,j-1)==2)) == size(ONorOFF,1)) % This checks to see if there is already a diode state for the next interval that does not affect switching actions
                                                [ ts, ~, ~, ~, ~,keep_SS] = obj.Baxter_adjustDiodeConduction(obj.Xs,j,i,max(obj.Xs(i,:)),-1,0.001,0,keep_SS);
                                                obj.setts(ts);
                                                order = obj.order;
                                                not_physical = true;
                                                obj.SS_Soln(1);
                                                obj.CorrectXs(1);
                                                
                                            else
                                                %                                                         replace = 0; new_index = 1;
                                                %                                                         [ts,order,new_index] = obj.adjust_time_single(sign,new_index,ts,order,waveform,i,k,time_ratio,ONorOFF,replace);
                                                %                                                         ts = [ts(1:k-1) ts(k)*.95 ts(k)*.05 ts(k+1:end)];
                                                %                                                         obj.setts(ts);
                                                %                                                         obj.SS_Soln(keep_SS);
                                                %                                                         obj.CorrectXs(keep_SS);
                                                %                                                         [ ts, ~, ~, ~, ~,keep_SS] = obj.Baxter_adjustDiodeConduction(obj.Xs,j,i,max(y(StateNumbers(i),:)),-1,0.001,0,keep_SS);
                                                %                                                         obj.setts(ts);
                                                %                                                         order = obj.order;
                                                %                                                         not_physical = true;
                                                %                                                         break
                                                
                                                last_violations_bd_turn_off(i,j-1) = 1;
                                                
                                            end
                                        end
                                    end
                                end
                                
                                if obj.Xs(i,j-1)>-1 % if there is only a violation at the beginning of the time interval then look to set the end of the last time interval to be equal to the diode forward votlage (1)
                                    first_violations = 1;
                                    
                                    if j == 2
                                        j = size(Xs,2);
                                        if first_violations %&& %ONorOFF(i,end) == 2
                                            
                                            %                                                     [ ts, ~, ~, ~, ~,keep_SS] = obj.Baxter_adjustDiodeConduction(obj.Xs,j-1,i,-1,min(obj.Xs(i,:)),0.001,1,keep_SS);
                                            order = obj.order;
                                            %                                                     not_physical = true;
                                            %                                                  %   obj.SS_Soln(1);
                                            %          obj.CorrectXs(1);
                                        end
                                        
                                    else
                                        
                                        if first_violations %&& %ONorOFF(i,k-1) == 2
                                            
                                            
                                            if ONorOFF(i,j-2) == -1 && (sum((ONorOFF(:,j-1)==2)==(ONorOFF(:,j-2)==2)) == size(ONorOFF,1)) % This checks to see if there is already a diode state for the next interval that does not affect switching actions
                                                [ ts, ~, ~, ~, ~,keep_SS] = obj.Baxter_adjustDiodeConduction(obj.Xs,j-1,i,-1,min(obj.Xs(i,:)),0.001,1,keep_SS);
                                                obj.setts(ts);
                                                order = obj.order;
                                                not_physical = true;
                                                %  obj.SS_Soln(1);
                                                %      obj.CorrectXs(1);
                                                
                                            else
                                                %                                                         replace = 0; new_index = 0;
                                                %                                                         [ts,order,new_index] = obj.adjust_time_single(sign,new_index,ts,order,waveform,i,k,time_ratio,ONorOFF,replace);
                                                %                                                         ts = [ts(1:k-1) ts(k)*.05 ts(k)*.95 ts(k+1:end)];
                                                %                                                         obj.setts(ts);
                                                %                                                         obj.SS_Soln(keep_SS);
                                                %                                                         obj.CorrectXs(keep_SS);
                                                %                                                         [ ts, ~, ~, ~, ~,keep_SS] = obj.Baxter_adjustDiodeConduction(obj.Xs,j,i,max(y(StateNumbers(i),:)),-1,0.001,0,keep_SS);
                                                %                                                         obj.setts(ts);
                                                %                                                         order = obj.order;
                                                %                                                         not_physical = true;
                                                %                                                         break
                                                
                                                first_violations_bd_turn_off(i,j-1) = 1;
                                                
                                            end
                                            
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                                
                                
                                %                                     else
                                %
                                %
                                %
                                %
                                %
                                %                                         [ts,order,new_index] = obj.adjust_time(sign,new_index,ts,order,waveform,i,k,time_ratio);
                                %
                                %                                         not_physical = true;
                                %                                     end
                                
                                
                                
                                
                                
                                %                                 if sum(waveform>-1+1*tol)>0 && debug
                                %                                     fprintf('Body Diode conducting: %s in time interval %.0f \n',obj.Converter.Topology.Parser.StateNames{i,1},j-1)
                                %                                     sign = [1,1];
                                %                                     [ts,order,new_index] = obj.adjust_time(sign,new_index,ts,order,waveform,i,k,time_ratio);
                                %                                     not_physical = true;
                                %                                     if waveform(1)>-1
                                %                                         first_violations = 1;
                                %                                     end
                                %                                     break
                                %                                 end
                                
                            end
                            
                        else
                            fprintf('Messed up\n')
                        end
                    else
                        % fprintf('Found a FET or Diode that wasn''t a FET or Diode\n Don''t panic\n')
                    end
                    
                end
                
            end
            if not_physical
                Jello = 8989;
            end
            j = j+1;
            
        end
        %             end
        if the_big_counter<2
            [ts,ONorOFF]=obj.adjust_time_single_all(ts,last_violations_bd_turn_on,last_violations_bd_turn_off, first_violations_bd_turn_on,  first_violations_bd_turn_off) ;
        end
        
        
        
        
        if not_physical % If there were no changes made then skip this step
            % Combine similar adjacent states would be nice!
            
            % This fixes the time intervals so there are not two adjacent time
            % intervals with the same state
            %                 the_size = length(order);
            %                 the_key = 1;
            %                 while the_key < the_size
            %                     if order(the_key)==order(the_key+1)
            %                         order(the_key+1) = [];
            %                         ts(the_key) = ts(the_key) + ts(the_key+1);
            %                         ts(the_key+1) = [];
            %                         the_size = the_size-1;
            %                         the_key = the_key-1;
            %                     end
            %                     the_key = the_key +1;
            %                 end
            %
            %
            %
            %                 obj.ts = ts;
            %                 obj.order = order;
            
            ONorOFF = obj.Converter.Topology.Parser.ONorOFF;
            time_intervals = size(ONorOFF,2);
            key = 1;
            while key < time_intervals
                
                if sum(ONorOFF(:,key)==ONorOFF(:,key+1))==size(ONorOFF,1) || ts(key+1) < 1e-12
                    obj.As(:,:,key+1) = [];
                    obj.Bs(:,:,key+1) = [];
                    obj.Cs(:,:,key+1) = [];
                    obj.Ds(:,:,key+1) = [];
                    ONorOFF(:,key+1) = [];
                    ts(key) = ts(key) + ts(key+1);
                    ts(key+1) = [];
                    
                    time_intervals = time_intervals-1;
                    key = key-1;
                    if debug
                        fprintf('Found a repeated state \n')
                    end
                end
                
                key = key+1;
            end
            obj.Converter.Topology.Parser.ONorOFF = ONorOFF;
            obj.setts(ts);
            
            
            %                 obj.updateTestConverter();
            obj.SS_Soln(keep_SS);
            obj.CorrectXs(keep_SS);
            % obj.Converter.Topology.Parser.find_diode_new(obj.order);
        end
        
        history_i(end+1)=i;
        history_j(end+1)=j;
        
        
    end
    check = 1;
    %obj.Xs(:,1) = obj.Xs(:,end); % Test step twords ss soln
    obj.SS_Soln();
    obj.CorrectXs();
    %obj.Converter.Topology.Parser.find_diode_new(obj.order);
    
    not_reached_SS = ~isequal(goal_SS,obj.Xs(:,1));
    %obj.Baxter_adjustDiodeConduction(Xs,3,2,1,-100)
    
    ts_og = ts;
    % 999 by 999 matrix takes over an hour to compute
    found_error = zeros(99,999);
    for primary = 1:1:99
        for secondary  = 1:1:99
            iter_pri = ts_og(2);
            iter_sec = ts_og(5);
            obj.setts([ts_og(1)-iter_pri*(primary-1) ts_og(2)+iter_pri*(primary-1) ts_og(3) ts_og(4)-iter_sec*(secondary-1) ts_og(5)+iter_sec*(secondary-1) ts_og(6) ts_og(1)-iter_pri*(primary-1) ts_og(2)+iter_pri*(primary-1) ts_og(3) ts_og(4)-iter_sec*(secondary-1) ts_og(5)+iter_sec*(secondary-1) ts_og(6)])
            obj.SS_Soln();
            obj.CorrectXs();
            primary_error = obj.Xs(2,2)+1;
            secondary_error = obj.Xs(4,5)+1;
            found_error(primary,secondary) =  norm([primary_error,secondary_error]);
        end
    end
    
        
figure(3)
surf([1:secondary]*iter_sec,[1:primary]*iter_pri,found_error)
        
        
   Alli = 124512315;    
    
    
    
    
catch ME
    rethrow(ME)
end



end % That's all Folks



