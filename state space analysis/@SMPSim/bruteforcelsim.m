function [check] = bruteforcelsim(obj,iterations)
%BRUTEFORCELSIM Attempts to optimize the time intervals of a converter using lsim()
%   iterations is the number of iterations that will be looped through to try and optimize the converter provided in the class
%
%     %%%%%%   %      %  %%%%%%%    %%%%%%
%    %      %  %      %  %      %  %      %
%    %      %  %      %  %      %  %      %
%    %%%%%%%%  %      %  %%%%%%%   %%%%%%%%
%    %      %  %      %  %%        %      %
%    %      %  %      %  % %       %      %
%    %      %  %      %  %  %      %      %
%    %      %  %      %  %   %     %      %
%    %      %   %    %   %    %    %      %
%    %      %    %%%%    %     %   %      %



% VDmax is the limited forward voltage allowed on a diode
% IDmax is the reverse current allowed on a diode
% VMmax is the limited forward voltage allowed on a diode
% IMmax is the reverse current allowed on a body diode (FET OFF)
VDmax = 3;
IDmax = -0.1;
VMmax = -3;
IMmax = 0.1;

i = [];

the_counter = 0;
Xs = obj.Xs;
while the_counter<=iterations
    
    the_counter = the_counter+1;
    
    % Important set up stuff
    debug = true;
    [xs, t, y, time_interval] = obj.SS_WF_Reconstruct();
    StateNumbers = obj.Converter.Topology.Parser.StateNumbers;
    StateNumbers_Opp = obj.Converter.Topology.Parser.StateNumbers_Opposite;
    ONorOFF = obj.Converter.Topology.Parser.ONorOFF;
    
    if debug
        fprintf('--------- \n')
        fprintf('Iteration number %.0f \n',the_counter)
        % What we are given
        figure
        ns = size(xs,1);
        StateNumbers = obj.Converter.Topology.Parser.StateNumbers;
        for z=1:ns
            subplot(10*ns,1,z*10-9:z*10)
            hold on;
            plot(t,y(StateNumbers(z),:), 'Linewidth', 3);
            ylabel(obj.getstatenames{z})
            box on
            
            if(z<ns)
                set(gca, 'Xticklabel', []);
            else
                xlabel('t(s)')
            end
        end
    end
    
    
    
    for i = 1:1:size(Xs,1) % Cycle through state variables
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
        
        %% Cycle through time intervals
        for j = 2:1:size(Xs,2)
            j;% is time interval for Xss
            k = j-1; % k is time interval for everything else
            if ONorOFF(i,k) ~=0 % if FET or Diode
                if obj.Converter.Topology.Parser.DMpos(i,2)==1 % if diode
                    % Determine if a state changes from the Diode
                    % being on to being off or vice vera.
                    if ONorOFF(i,k) == 1 % if diode ON
                        if ~isempty(find(y(StateNumbers(i),time_interval(i)+1:time_interval(i+1)) < 0,1)) && debug
                            fprintf('State Violation (Diode turn off) of %s in time interval %.0f \n',obj.Converter.Topology.Parser.StateNames{i,1},j-1)
                        end
                    elseif ONorOFF(i,j-1) == -1 % if diode off
                        if ~isempty(find(y(StateNumbers_Opp(i),time_interval(i)+1:time_interval(i+1)) > 0,1)) && debug
                            fprintf('State Violation (Diode turn on) of %s in time interval %.0f \n',obj.Converter.Topology.Parser.StateNames{i,1},j-1)
                        end
                    else
                        fprintf('Messed up\n')
                    end
                elseif obj.Converter.Topology.Parser.DMpos(i,3)==1 % if FET
                    if ONorOFF(i,j-1) == 1 % if FET ON
                        if ~isempty(find(y(StateNumbers_Opp(i),time_interval(i)+1:time_interval(i+1)) < 0,1)) && debug
                            % Check if body diode conducts (time interval could be shortened)
                            fprintf('Body Diode conducting: %s in time interval %.0f \n',obj.Converter.Topology.Parser.StateNames{i,1},j-1)
                        end
                        
                    elseif ONorOFF(i,j-1) == -1 % if FET off
                        % Can not have positive ID
                        if ~isempty(find(y(StateNumbers_Opp(i),time_interval(i)+1:time_interval(i+1)) > IMmax,1)) && debug % Only finds first violation to improve speed
                            fprintf('State Violation of FET reverse current %s exceed %.2f A \n',obj.Converter.Topology.Parser.StateNames{i,1},IDmax)
                        end
                        
                        if ~isempty(find(y(StateNumbers_Opp(i),time_interval(i)+1:time_interval(i+1)) < 0,1)) && debug
                            fprintf('Body Diode conducting: %s in time interval %.0f \n',obj.Converter.Topology.Parser.StateNames{i,1},j-1)
                        end
                        
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
    end
    
    
    
    % Determine if it its possible to reach ideal value for soft switching (lsim for length of period)
    % If it is possible to reach value then set the time to the first instance of that value
    % If is not possible to reach value then step in direction by a certain percentage of initial guess.
    if mod(the_counter,2)
        %% Adjust Deadtimes
        deadtimes = obj.dead_time_intervals ;
        deadstates = obj.dead_time_states;
        deadgoals = obj.dead_time_goals;
        ts = obj.ts;
        Ts = sum(ts);
        time_ratio = Ts/time_interval(end);
        
        % How do I translate between time in the SS_WF domain to the acutal
        % time lengths of each state that I need to set to run SSsoln.
        
        for i = 1:1:length(deadtimes)
            index = time_interval(deadtimes(i))+1:time_interval(deadtimes(i)+1);
            waveform = round(y(StateNumbers(deadstates(i)),index),0);
            [V,P] = find(waveform == deadgoals(i));
            
            % If it is not already fairly close
            if ~isempty(V)
                ts(deadtimes(i))=time_ratio*P(1);
                if size(obj.As,3)==deadtimes(i)
                    ts(1) = ts(1)+time_ratio*(time_interval(deadtimes(i)+1)-time_interval(deadtimes(i))-P(1));
                else
                    ts(deadtimes(i)+1) = ts(deadtimes(i)+1)+time_ratio*(time_interval(deadtimes(i)+1)-time_interval(deadtimes(i))-P(1));
                end
                % Set time for this interval to be this value if overshoot
            else
                % If not overshoot then test
                
                [V,P]=min(abs(deadgoals(i)-round(waveform,3)));
                ts(deadtimes(i))=time_ratio*P(1);
                
                if size(obj.As,3)==deadtimes(i)
                    ts(1) = ts(1)+time_ratio*(time_interval(deadtimes(i)+1)-time_interval(deadtimes(i))-P(1));
                else
                    ts(deadtimes(i)+1) = ts(deadtimes(i)+1)+time_ratio*(time_interval(deadtimes(i)+1)-time_interval(deadtimes(i))-P(1));
                end
                % Move towards the region of less error (reduce error)
            end
            
        end
   
        
        
        
    else
         %% Adjust Power times
        
        Vo_index = obj.Vo_index;
        Vo_ideal_value = obj.Vo_ideal_value;
        average = mean(y(StateNumbers(Vo_index),:));
        Voerr = Vo_ideal_value-average;
        Perturb1_index = obj.Perturb1_index;
        Perturb2_index = obj.Perturb2_index;
        
        
        delta_DTs = max(min(ts)/10, sum(ts)/10000);
        [dXs]=obj.Baxter_StateSensitivity('ts', Perturb1_index, delta_DTs,Perturb2_index);
        dxsdt = (dXs-obj.Xs)/delta_DTs;
        dt = (Voerr)/mean(dxsdt(Vo_index,:));
        
        ts(1) = ts(1)+0.5*dt;
        ts(3) = ts(3)-0.5*dt;
        
    end
    
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
    
    obj.ts = ts;
    obj.SS_Soln();
    obj.CorrectXs();
    
    
end
check =1;
end % That's all Folks
