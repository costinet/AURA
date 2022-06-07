function [not_reached_SS] = Three_tier_diode_correct_num(obj,iterations,debug,re_start)
%THREE_TIER_Jacobian is a test of the new jacobian method of solving for steady state
%takes the steady state solution of a
%converter can checks to determine if there are any diode violations
%and correct them

%     _   _   _  ____    _
%    / \ | | | |/ _  |  / \
%   / _ \| | | | (_| | / _ \
%  / ___ | |_| |> _  |/ ___ \
% /_/   \_\___//_/ |_/_/   \_\

% Set debug flag

%
%obj.Xs(:,end) = [0;5;5;0;0;0;0];
% obj.Xs(:,end) = [0;0;0;0;0;0;0];
% obj.Xs(:,end) = [0;0;0;0;0];


if re_start
    if ~isempty(obj.As_saved)
        obj.As = obj.As_saved ;
        obj.Bs = obj.Bs_saved ;
        obj.Cs = obj.Cs_saved ;
        obj.Ds = obj.Ds_saved ;
        obj.u = obj.u_saved;
        obj.eigA = obj.eigA_saved;
        obj.Converter.Topology.Parser.ONorOFF = obj.ONorOFF_saved;
        obj.ts = obj.ts_saved;
        obj.Xs = obj.Xs_saved;
    end
    
    
else
    % Set this for the lsim function to be able to run without issue
    obj.As_OG = obj.As;
    obj.Bs_OG = obj.Bs;
    obj.Cs_OG = obj.Cs;
    obj.Ds_OG = obj.Ds;
    obj.u_OG = obj.u;
    obj.eigA_OG = obj.eigA;
    obj.ONorOFF_OG = obj.Converter.Topology.Parser.ONorOFF;
    obj.ts_OG = obj.ts;
end






gdiode = zeros(size(obj.Converter.Topology.Parser.ONorOFF));
ts_hist = [];
Xs_hist = [];
Time_stamp = [];
% Initalization of Variables and Counters
the_big_counter = 0;
more_iterations=iterations;
not_reached_SS = true;
previous_multi_violation = false;
breakbreak = false;
tol = 0.1;
keep_SS = false; % dont really use this anymore but it is still riddled throughout the code

Vf = obj.Converter.Topology.Parser.Fwd_Voltage;
bad_converter = 0;


if debug
    %Plot_Waveforms;
end

while not_reached_SS && the_big_counter<=more_iterations
    
    if the_big_counter >15
        that_a_lot = 5645648456;
    end
    
    if re_start && the_big_counter ==2
        
        obj.As_saved = obj.As  ;
        obj.Bs_saved = obj.Bs ;
        obj.Cs_saved = obj.Cs ;
        obj.Ds_saved = obj.Ds ;
        obj.u_saved = obj.u;
        obj.eigA_saved = obj.eigA;
        obj.ONorOFF_saved = obj.Converter.Topology.Parser.ONorOFF ;
        obj.ts_saved = obj.ts;
        obj.Xs_saved = obj.Xs;
        
    end
    
    
    
    the_big_counter = the_big_counter+1;
    
    try % For debugging
        not_reached_SS = false;
        %% Eigenvalue check
        
        % The eigenvalue check is done at the beginning of the function to
        % determine which
        
        
        
        
        imagine_eigA = imag(obj.eigA); % get the imaginary part of the eigenvalue
        imagine_eigA(imagine_eigA(:)<=0) = 0; % Only take the postive part of the eigenvalue
        
        num_eigA_volations=sum(1./obj.ts<imagine_eigA./pi(),1); % Find the time intervals that have an eigenvalue that has faster dynamics than its length
        
        multi_violations = sum(num_eigA_volations>1)>0;
        
        %    if multi_violations == 1
        %        bad_converter = 1;
        %    end
        
        %   if the_big_counter == 10 || the_big_counter == 11 || the_big_counter == 20 || the_big_counter == 21 || the_big_counter == 30 || the_big_counter == 31
        %       multi_violations = 1;
        %   end
        
        
        if re_start
            %  For flyback
           %  if the_big_counter < 2 && mod(the_big_counter,2)
            %     multi_violations = 1;
            % end
        end
        
        if multi_violations
            % Kick to code that will lsim throught the code if there is more
            % than one eigenvalue that is faster than its corresponding time interval
            [~,obj.Xs] = obj.PLECS_lsim_solve(obj.Xs(:,end),1);
            gdiode=zeros(size(obj.Converter.Topology.Parser.ONorOFF));
            obj.ts_history = [];
            obj.Xs_history = [];
            previous_multi_violation = true;
            not_reached_SS = true;
            disp('Multi_Violations time step solve');
            continue
        elseif previous_multi_violation
            obj.SS_Soln();
            obj.CorrectXs();
            previous_multi_violation = false;
            not_reached_SS = true;
            continue
        end
        
        
        
        
        
        if debug
            %Plot_Waveforms;
        end
        
        %% State Sensitvity Check
        
        
        Xs = obj.Xs;
        ONorOFF = obj.Converter.Topology.Parser.ONorOFF;
        ts = obj.ts;
        
        last_violations_bd_turn_on = zeros(size(ONorOFF));
        last_violations_bd_turn_off = last_violations_bd_turn_on;
        first_violations_bd_turn_off = last_violations_bd_turn_off;
        first_violations_bd_turn_on = last_violations_bd_turn_off;
        
        
        j = 2;
        time_variable_size = size(Xs,2);
        
        ONorOFF_ON = ONorOFF == 2;
        ONorOFF_OFF = ONorOFF == -1;
        ONorOFF_D = ONorOFF == 1;
        
        Vf_compare = repmat(-Vf,1,size(Xs,2));
        
        Violation_B_D = Xs(:,1:end-1)>Vf_compare(:,1:end-1)&ONorOFF_D;
        Violation_B_OFF = Xs(:,1:end-1)<Vf_compare(:,1:end-1)&ONorOFF_OFF;
        
        Violation_A_D = Xs(:,2:end)>Vf_compare(:,2:end)&ONorOFF_D;
        Violation_A_OFF = Xs(:,2:end)<Vf_compare(:,2:end)&ONorOFF_OFF;
        
        
        
        if sum(Violation_B_D + Violation_B_OFF+ Violation_A_D +Violation_A_OFF,'all')==0
            % there is no violation
            return
        else
            J = 456456465;
        end
        altered = false;
        last_violations_bd_turn_on = circshift(ONorOFF,-1,2) .* Violation_A_OFF == -1;
        if sum(last_violations_bd_turn_on,'all')>0
            altered = true;
        [ts,ONorOFF,gdiode]=obj.adjust_time_single_all_num(ts,last_violations_bd_turn_on,last_violations_bd_turn_off, first_violations_bd_turn_on,  first_violations_bd_turn_off,gdiode) ;
        not_reached_SS = true;
        end
        
        if ~altered
            
            
          
            
            
            order  = 1;
            [J1, J2, XssF, XssB, X0, dt] = obj.Baxter_Jacobianfunction(order);
            
            
            (Violation_A_OFF.*(Vf_compare(:,2:end)-Xs(:,2:end)))./J1
            
            
        end
        
        %if not_physical % If there were no changes made then skip this step
        % Combine similar adjacent states would be nice!ex
        
        ONorOFF = obj.Converter.Topology.Parser.ONorOFF;
        time_intervals = size(ONorOFF,2);
        key = 1;
        while key < time_intervals
            
            if sum(ONorOFF(:,key)==ONorOFF(:,key+1))==size(ONorOFF,1) || ts(key+1) < 1e-12
                obj.As(:,:,key+1) = [];
                obj.Bs(:,:,key+1) = [];
                obj.Cs(:,:,key+1) = [];
                obj.Ds(:,:,key+1) = [];
                obj.eigA(:,key+1) = [];
                ONorOFF(:,key+1) = [];
                gdiode(:,key+1) = [];
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
        
        %end
        
        
        breakbreak = 0;
        obj.SS_Soln();
        obj.CorrectXs();
        
        obj.ts_hist(obj.ts);
        obj.Xs_hist(obj.Xs);
        %    Time_stamp(end+1) = toc;
    catch ME
        ME
        ME.stack.line
        rethrow(ME)
    end
end

%if not_reached_SS
%    fprintf('Failed to converge to steady state solution \n')
%end



%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Space to work on state sensitivity issues
vector = linspace(1e-12,2.9999e-08,1000);

%%{
for i = 1:length(vector)
    obj.ts(11) = vector(i);
    obj.ts(12) = 3.0000e-08-vector(i);
    
    obj.SS_Soln();
    obj.CorrectXs();
    Xs_history_example(:,:,i) = obj.Xs;
    % end
    
    %{
    figure
    ns = 7;
    for z=1:ns
        ax = subplot(10*ns,1,z*10-9:z*10);
        hold on;
        plot(vector,squeeze(Xs_history_example(z,6,:)), 'Linewidth', 3);
        ylabel(obj.getstatenames{z})
        box on
        if(z<ns)
            set(gca, 'Xticklabel', []);
        else
            xlabel('time of Ti=11')
        end
    end
    
    %}
    
    delta_DTs = 10e-12;  % 10e-14;
    keep_SS = false;
    
    %vector = linspace(delta_DTs,delta_DTs*1000,1000);
    %for i = 1:length(vector)
    [dXs,delta_DTs] = obj.Baxter_StateSensitivity2(keep_SS, 'ts', 11, delta_DTs, 12);
    dxsdt(:,:,i) = (dXs-obj.Xs)/delta_DTs;
    
    [dXs_back,delta_DTs_back] = obj.Baxter_StateSensitivity2(keep_SS, 'ts', 12, delta_DTs, 11);
    dxsdt_back(:,:,i) = (obj.Xs - dXs_back)/delta_DTs_back;
    
    Second_Derivative(:,:,i) = (dXs_back-2*obj.Xs+dXs)/(delta_DTs^2);
    
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





    function Plot_Waveforms(statenum,oppstatenum,plotxparam)
        %PLOT_WAVEFROMS is a nested function that plots all of the states
        %and the inverse states (V->I or I->V). They will appear in
        %figures 100 and 101
        
        
        [xs, t, y, time_interval] = obj.SS_WF_Reconstruct();
        StateNumbers = obj.Converter.Topology.Parser.StateNumbers;
        StateNumbers_Opp = obj.Converter.Topology.Parser.StateNumbers_Opposite;
        
        switch nargin
            case 0
                statenum = 100;
                oppstatenum = 101;
                plotxparam  = 0;
            case 1
                oppstatenum = 101;
                plotxparam = 0;
            case 2
                plotxparam = plotxparam;
        end
        
        if plotxparam~=0
            ns = length(plotxparam);
            figure(50)
            for z=1:ns
                ax = subplot(10*ns,1,z*10-9:z*10);
                
                hold on;
                if plotxparam(z)<0
                    plotxparam(z) = abs(plotxparam(z));
                    plot(t*10^6,-y(plotxparam(z),:), 'Linewidth', 3);
                else
                    plot(t*10^6,y(plotxparam(z),:), 'Linewidth', 3);
                end
                ylabel(plotxparam(z))
                box on
                %ax.YLim = [min(xs(z,:))-abs(0.5*min(xs(z,:))) max(xs(z,:))+abs(0.5*max(xs(z,:)))];
                if(z<ns)
                    set(gca, 'Xticklabel', []);
                else
                    xlabel('t [\mus]')
                end
            end
        end
        
        
        figure(statenum)
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
        
        
        figure(oppstatenum)
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
% Save function would go here
end