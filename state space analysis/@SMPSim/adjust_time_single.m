function [ts,order,new_index] = adjust_time_single(obj,sign,new_index,ts,order,waveform,i,k,time_ratio,ONorOFF,replace)
%adjaust_time calcualtes the time needed and ordered states needed to
%adjust to have a correct physical circuit with diode implementation
%   Detailed explanation goes here


% sign is a 1x2 double that determines the states that are chosen to
% replace part or all of the period chosen. sign(1) determines whether
% there is a < or > sign used to determine violations (could probably
% be reduced with *(-1)). sign(2) determines the sign of the sign of
% the forward voltage to make the comparison. This can be changed once
% we get in a database of component values in a class.


%     _   _   _  ____    _
%    / \ | | | |/ _  |  / \
%   / _ \| | | | (_| | / _ \
%  / ___ | |_| |> _  |/ ___ \
% /_/   \_\___//_/ |_/_/   \_\


% Add one variable at the end of a time interval

try
    % OnorOFF is organized so states are the rows and time intervals
    % are the columns
    %     binary = obj.binary;
    
    %     binary(binary==1)=2; % Turn logic 1 (ON) to a 2
    %     binary(binary==0)=-1; % Turn logic 0 (OFF) to a -1
    
    if replace
        % This is for when the error is for the entire time interval
        % so there is no additional time interval added to the system,
        % instead the previous time interval is adjusted to account
        % for the error 
        new_state = ONorOFF(:,k);
        new_state(i,1) = -1*new_state(i,1);
        new_state_space = [ONorOFF(:,1:k-1), new_state ,ONorOFF(:,k+1:end) ];
        obj.Converter.Topology.Parser.ONorOFF = new_state_space; % Adjustment of ONorOFF matrix
        [A,B,C,D] = obj.add_state_matrix(new_state_space,k);
        
        new_A(:,:,1:k-1) = obj.As(:,:,1:k-1);
        new_A(:,:,k) =  A;
        new_A(:,:,end+1:size(obj.As,3)) = obj.As(:,:,k+1:end) ;
        new_B(:,:,1:k-1) = obj.Bs(:,:,1:k-1);
        new_B(:,:,k) =  B;
        new_B(:,:,end+1:size(obj.Bs,3)) = obj.Bs(:,:,k+1:end) ;
        new_C(:,:,1:k-1) = obj.Cs(:,:,1:k-1);
        new_C(:,:,k) =  C;
        new_C(:,:,end+1:size(obj.Cs,3)) = obj.Cs(:,:,k+1:end) ;
        new_D(:,:,1:k-1) = obj.Ds(:,:,1:k-1);
        new_D(:,:,k) =  D;
        new_D(:,:,end+1:size(obj.Ds,3)) = obj.Ds(:,:,k+1:end) ;
        
    else
        if new_index
            % This is for when the error occurs at the end of
            % the time interval and thus the new diode state should be
            % after the time interval where the error was observed
            new_state = ONorOFF(:,k);
            new_state(i,1) = -1*new_state(i,1);
            new_state_space = [ONorOFF(:,1:k), new_state ,ONorOFF(:,k+1:end) ];
            obj.Converter.Topology.Parser.ONorOFF = new_state_space;
            [A,B,C,D] = obj.add_state_matrix(new_state_space,k+1);
            
            new_A(:,:,1:k) = obj.As(:,:,1:k);
            new_A(:,:,k+1) =  A;
            new_A(:,:,end+1:size(obj.As,3)+1) = obj.As(:,:,k+1:end) ;
            new_B(:,:,1:k) = obj.Bs(:,:,1:k);
            new_B(:,:,k+1) =  B;
            new_B(:,:,end+1:size(obj.Bs,3)+1) = obj.Bs(:,:,k+1:end) ;
            new_C(:,:,1:k) = obj.Cs(:,:,1:k);
            new_C(:,:,k+1) =  C;
            new_C(:,:,end+1:size(obj.Cs,3)+1) = obj.Cs(:,:,k+1:end) ;
            new_D(:,:,1:k) = obj.Ds(:,:,1:k);
            new_D(:,:,k+1) =  D;
            new_D(:,:,end+1:size(obj.Ds,3)+1) = obj.Ds(:,:,k+1:end) ;
            
        else
            
            % This is for when the error occurs at the beginning of
            % the time interval and thus the new diode state should be
            % after the time interval where the error was observed
            
            new_state = ONorOFF(:,k);
            new_state(i,1) = -1*new_state(i,1);
            new_state_space = [ONorOFF(:,1:k-1), new_state ,ONorOFF(:,k:end) ];
            obj.Converter.Topology.Parser.ONorOFF = new_state_space;
            [A,B,C,D] = obj.add_state_matrix(new_state_space,k);
            
            new_A(:,:,1:k-1) = obj.As(:,:,1:k-1);
            new_A(:,:,k) =  A;
            new_A(:,:,end+1:size(obj.As,3)+1) = obj.As(:,:,k:end) ;
            new_B(:,:,1:k-1) = obj.Bs(:,:,1:k-1);
            new_B(:,:,k) =  B;
            new_B(:,:,end+1:size(obj.Bs,3)+1) = obj.Bs(:,:,k:end) ;
            new_C(:,:,1:k-1) = obj.Cs(:,:,1:k-1);
            new_C(:,:,k) =  C;
            new_C(:,:,end+1:size(obj.Cs,3)+1) = obj.Cs(:,:,k:end) ;
            new_D(:,:,1:k-1) = obj.Ds(:,:,1:k-1);
            new_D(:,:,k) =  D;
            new_D(:,:,end+1:size(obj.Ds,3)+1) = obj.Ds(:,:,k:end) ;
            
            
            
        end
        
    end
    obj.As = new_A;
    obj.Bs = new_B;
    obj.Cs = new_C;
    obj.Ds = new_D;
    
    
    
    
    % if sign(1)  == 1
    %     [P] = find(diff(waveform>1*(-1)^sign(2))~=0); % Find all cases of violations
    %     flipflop = waveform(1)>1*(-1)^sign(2); % Check to see if the inital value was in violation. This will be used later on if there are multiple crossings of violations within the same time interval. We simply flipflop between states in the order variable
    %     ts_new = zeros(1,length(P)+1); % intialize new ts two be inserted for the time
    %     order_new = ts_new;
    %     % Get all of the values we need from Parser
    %     if sign(2) == 1
    %         bd_state = obj.Converter.Topology.Parser.BD_OFF_state;
    %     end
    %     if sign(2) == 0
    %         bd_state = obj.Converter.Topology.Parser.BD_state;
    %     end
    %     old_state = obj.order(k);
    %     ordered_state_index = obj.Converter.Topology.Parser.OrderedNamesnum;
    %     switches = obj.Converter.Topology.Parser.Switches;
    %     bd_state_new = bd_state(obj.order(k),ordered_state_index(i,obj.order(k))==switches);
    %
    %     if isempty(P)
    %         % if there is no crossing then the entire state is in
    %         % violation of the
    %         order(new_index) = bd_state_new(1);
    %     else
    %         for number_of_changes = 1:1:length(P)+1
    %             order_new(number_of_changes) = flipflop;
    %             if number_of_changes == 1
    %                 ts_new(1) = time_ratio*(P(1));
    %
    %             elseif number_of_changes == length(P)+1
    %                 ts_new(number_of_changes) = time_ratio*(length(waveform)-P(number_of_changes-1));
    %
    %             else
    %                 ts_new(number_of_changes) = time_ratio*(P(number_of_changes)-P(number_of_changes-1));
    %
    %             end
    %
    %             flipflop=~flipflop;
    %         end
    %
    %         order_new(order_new==1) = bd_state_new(1);
    %         order_new(order_new==0) = old_state;
    %         ts = [ts(1:new_index-1) ts_new ts(new_index+1:end)];
    %         order = [order(1:new_index-1) order_new order(new_index+1:end)]; % Assume that FET is off and body diode is ON will have to adjust then when that correction is made its the (1) index
    %
    %         new_index = new_index+length(P);
    %     end
    % end
    %
    % if sign(1)  == 0
    %     [P] = find(diff(waveform<1*(-1)^sign(2))~=0);
    %     flipflop = waveform(1)<1*(-1)^sign(2);
    %     ts_new = zeros(1,length(P)+1);
    %     order_new = ts_new;
    %
    %     if sign(2) == 0
    %         bd_state = obj.Converter.Topology.Parser.BD_OFF_state;
    %     end
    %     if sign(2) == 1
    %         bd_state = obj.Converter.Topology.Parser.BD_state;
    %     end
    %
    %     old_state = obj.order(k);
    %     ordered_state_index = obj.Converter.Topology.Parser.OrderedNamesnum;
    %     switches = obj.Converter.Topology.Parser.Switches;
    %     bd_state_new = bd_state(obj.order(k),ordered_state_index(i,obj.order(k))==switches);
    %     if isempty(P)
    %         order(new_index) = bd_state_new(1);
    %     else
    %         for number_of_changes = 1:1:length(P)+1
    %             order_new(number_of_changes) = flipflop;
    %             if number_of_changes == 1
    %                 ts_new(1) = time_ratio*(P(1));
    %
    %             elseif number_of_changes == length(P)+1
    %                 ts_new(number_of_changes) = time_ratio*(length(waveform)-P(number_of_changes-1));
    %
    %             else
    %                 ts_new(number_of_changes) = time_ratio*(P(number_of_changes)-P(number_of_changes-1));
    %
    %             end
    %
    %             flipflop=~flipflop;
    %         end
    %
    %         order_new(order_new==1) = bd_state_new(1);
    %         order_new(order_new==0) = old_state;
    %         ts = [ts(1:new_index-1) ts_new ts(new_index+1:end)];
    %         order = [order(1:new_index-1) order_new order(new_index+1:end)]; % Assume that FET is off and body diode is ON will have to adjust then when that correction is made its the (1) index
    %
    %         new_index = new_index+length(P);
    %
    %     end
    % end
    %
    %
    %
    
    
    
    
    
    
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
    
catch ME
    
    rethrow(ME)
end

end % That's all Folks

