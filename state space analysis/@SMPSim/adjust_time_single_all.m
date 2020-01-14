function [ts,new_state_space] = adjust_time_single_all(obj,ts,L_ON,L_OFF,F_ON,F_OFF)
%adjaust_time calcualtes the time needed and ordered states needed to
%adjust to have a correct physical circuit with diode implementation
%   Detailed explanation goes here


% Need design a sorting method to determine how to sort and create a
% be


%     _   _   _  ____    _
%    / \ | | | |/ _  |  / \
%   / _ \| | | | (_| | / _ \
%  / ___ | |_| |> _  |/ ___ \
% /_/   \_\___//_/ |_/_/   \_\


try
    
    % Add one variable at the end of a time interval
    ONorOFF = obj.Converter.Topology.Parser.ONorOFF;
    new_state_space = ONorOFF;
    tot_L_ON = sum(L_ON,1);
    tot_L_OFF = sum(L_OFF,1);
    tot_F_ON = sum(F_ON,1);
    tot_F_OFF = sum(F_OFF,1);
    
    new_A = obj.As;
    new_B = obj.Bs;
    new_C = obj.Cs;
    new_D = obj.Ds;
    new_eigA = obj.eigA;
    
    count = 1;
    k = 1;
    while k <=size(ONorOFF,2)
        
        if tot_L_ON(k)~=0
            
            new_state = ONorOFF(:,k); % Find state that needs to be copied
            new_state(L_ON(:,k)==1) = 1; % make correction in state
            
            Ti = k; % Current index of time interval
            Tc = k+1; % Next time interval
            
            % If the next time interval reaches the end of the defined period
            % then restart at beginning (1)
            if(Tc == size(ONorOFF,2) + 1)
                Tc = 1;
            end
            
            
            if new_state == ONorOFF(:,Tc)
                
                % This is where we can make a correction to determine
                % if there needs to be a separation of violations
                
                % Perform state sensitivity to see if there are values
                % of the violation that are common to a purtibation
                
                delta_DTs = max(min(ts)/1000, sum(ts)/100000);
                [dXs,delta_DTs] = obj.Baxter_StateSensitivity2(0, 'ts', Ti, delta_DTs, Tc);
                
                dxsdt = (dXs-obj.Xs)/delta_DTs;
                
                dxsdt_col = dxsdt(:,Tc).*L_ON(:,Ti);
                
                
                % use Repmat to create a column and row vector of the
                % values at the point of diode violation
                
                % Am looking to get something like a jacombian matrix
                % but its will be a set of differences
                
                jac_like = abs((repmat(dxsdt_col,[1,size(dxsdt,1)])-repmat(dxsdt_col',[size(dxsdt,1),1]))./repmat(dxsdt_col,[1,size(dxsdt,1)]));
                
                [Xcord,Ycord]=find(jac_like<0.01& jac_like~=0);
                
                
                jac_like = zeros(size(jac_like));
                All_cord=sub2ind(size(jac_like),Xcord,Ycord);
                jac_like(All_cord) = 1;
                
                tdelta = zeros(size(L_ON,1),1);
                % If this is true then there is a
                if true == all(all(jac_like==jac_like'))
                    for index = 1:1:size(L_ON,1)
                        
                        if L_ON(index,k)==1
                            tdelta(index) = (-1-obj.Xs(index,k+1))/(dxsdt(index,k+1));
                        end
                    end
                    [~,P]=min(tdelta);
                    L_ON(:,k) = jac_like(:,P);
                    L_ON(P,k) = 1;
                    k = k-1;
                    
                else
                    
                    fprintf('Not sure what to do here yet');
                end

                %{
                % match is of the form: columns correspond to Xcord and
                % rows correspond to Ycord
                for i = 1:1:length(Xcord)
                    match(:,i)=(Xcord==repmat(Ycord(i),[length(Xcord),1]))&(Ycord==repmat(Xcord(i),[length(Ycord),1]));
                    
                end
                 
                % If the transpose of the values are == 1 then both
                % are within the bounds and those two values can be
                % linked
                
                % If every element that was a match with its
                % reciprocole then we can shorthand get all of the
                % matches by just looking at the row and columns
                % quickly
                
                if true == all(all(match==match'))
                    index = 1;
                    together = zeros(size(Xcord,1));
                    while index < size(match,2)
                        match(:,match(index,:))=0;
                        
                        
                    end
                    
                end
               
                %}
                count = count-1;
            else
                new_state_space = [new_state_space(:,1:count), new_state ,new_state_space(:,count+1:end) ]; % Place new state space in the ONorOFF matrix
                [A,B,C,D,eigA] = obj.add_state_matrix(new_state); % Calculate the needed
                
                
                % ts = [ts(1:count-1) ts(count)*.05 ts(count)*.95 ts(count+1:end)];
                % This is the new one to find a finner mesh grid on the
                % surface plot
                ts = [ts(1:count-1) ts(count)*.99 ts(count)*.01 ts(count+1:end)];
                
                % This is the normal one:
                % ts = [ts(1:count-1) ts(count)*1.9074e-07/ts(count) ts(count)*(ts(count)-1.9074e-07)/ts(count) ts(count+1:end)];
                
                new_A(:,:,1:count) = new_A(:,:,1:count);
                new_A(:,:,count+2:size(new_A,3)+1) = new_A(:,:,count+1:end);
                new_A(:,:,count+1) =  A;
                new_B(:,:,1:count) = new_B(:,:,1:count);
                new_B(:,:,count+2:size(new_B,3)+1) = new_B(:,:,count+1:end) ;
                new_B(:,:,count+1) =  B;
                new_C(:,:,1:count) = new_C(:,:,1:count);
                new_C(:,:,count+2:size(new_C,3)+1) = new_C(:,:,count+1:end) ;
                new_C(:,:,count+1) =  C;
                new_D(:,:,1:count) = new_D(:,:,1:count);
                new_D(:,:,count+2:size(new_D,3)+1) = new_D(:,:,count+1:end) ;
                new_D(:,:,count+1) =  D;
                new_eigA(:,1:count) = new_eigA(:,1:count);
                new_eigA(:,count+2:size(new_eigA,2)+1) = new_eigA(:,count+1:end) ;
                new_eigA(:,count+1) =  eigA;
                
                count = count+1;
            end
            
        end
        
        count = count+1;
        k = k+1
    end
    
    obj.Converter.Topology.Parser.ONorOFF = new_state_space;
    obj.setts(ts);
    obj.As = new_A;
    obj.Bs = new_B;
    obj.Cs = new_C;
    obj.Ds = new_D;
    obj.eigA = new_eigA;
    
    
    
    
    
    
    
    
    %{
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
    %}
catch ME
    
    rethrow(ME)
end

end % That's all Folks

