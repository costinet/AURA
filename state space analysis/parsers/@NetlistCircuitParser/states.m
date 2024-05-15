function [NL,NLnets,forward_pass] = states(obj,state,i,switches)
% STATES finds the different switch positions of the netlist file converter
% and sets the net list for nodeloop by changing open switches to caps and
% on switches to caps in parallel with an on resistance.


%     _   _   _  ____    _    
%    / \ | | | |/ _  |  / \   
%   / _ \| | | | (_| | / _ \  
%  / ___ | |_| |> _  |/ ___ \ 
% /_/   \_\___//_/ |_/_/   \_\



% Does the current implementation have a forward voltage state and a
% non forward voltage state for diodes????

% Index values for components (will pass as variable... eventually):
numV = 1;
numBV = 2;
numMV = 3;
numC = 4;
numR = 5;
numL = 6;
numMI = 7;
numBI = 8;
numI = 9;
numD = 10;
numM = 11;

NL = obj.NL;
NLnets = obj.NLnets;
ON_States = cell(size(switches,2),1);
OFF_States = cell(size(switches,2),1);
ON_index = 1;
OFF_index = 1;
count  = 1;
forward_pass = zeros(size(obj.Diodes,1),1);
try
for j = 1:1:length(switches)
    
    % OPEN switch is replaced with cap
    if j~=length(switches)
        
        % if a FET then need to check that state and next state is
        % zero AND check that the current and next states correspond
        % to the same switch AND that is is a FET
        
        if (state(i,j) == 0 && state(i,j+1) == 0 && switches(j+1) == switches(j) && (NL(switches(j),1)==numM))
            row = zeros(size(NL,1),1); % Find element that we are trying to measure the voltage of
            row(switches(j)) = 1;
            
            
            % Add resistor in parallel to Cap:
            % Add new voltage measurement element to NLnets:
            NLnets(end+1,:) = [strcat(NLnets(switches(j),1),'_R_OFF'), NLnets(switches(j),2:end)];
            % Add new voltage measurement (current source) element to NL:
            NL(end+1,:)= [numR,sum(row.*NL(:,2:3)),NL(end,4)+1];
            
            % Set C to be current position:
            % Add new voltage measurement element to NLnets:
            NLnets(switches(j),:) = [strcat(NLnets(switches(j),1),'_C'), NLnets(switches(j),2:end)];
            % Add new voltage measurement (current source) element to NL:
            NL(switches(j),1)= [numC];
            
            % Add voltage source
            row = zeros(size(NL,1),1); % Find element that we are trying to measure the voltage of
            row(end) = 1;
            
            new_node = max(max(NL(:,2:3)))+1; % Set new node number (highest node number + 1)
            node = NL(sum(row.*NL(:,4)),2); % Find node number from element being measured
            
            New_Node = strcat(NLnets(sum(row.*NL(:,4)),2),'_Diode'); % Set new node cell name (highest node number + 1)
            Node = NLnets(sum(row.*NL(:,4)),2); % Find cell name from element being measured
            
            NLnets(sum(row.*NL(:,4)),2) = New_Node; % Adjust element beign measured node (cell net list)
            NLnets(end+1,:) = [strcat(NLnets(switches(j),1),'_VF'),New_Node,Node,NLnets(sum(row.*NL(:,4)),4:end)]; % Create measurement node (cell net list)
            %NLnets(sum(row.*NL(:,4)),1) = strcat(NLnets(switches(j),1),'_R_OFF');
            
            NL(sum(row.*NL(:,4)),2)=new_node; % Adjust element beign measured node (numerical net list)
            NL(end+1,:)= [numV,new_node,node,NL(end,4)+1]; % Create measurement node (numerical net list)
            ON_States(ON_index,1) = NLnets(end-1,1);
            ON_index = ON_index+1;
            count = count+1;
        end
        
        % CLOSED switch is replaced with cap and resistor
        if (state(i,j) == 1 && switches(j+1) == switches(j)) && (NL(switches(j),1)==numM)
            row = zeros(size(NL,1),1); % Find element that we are trying to measure the voltage of
            row(switches(j)) = 1;
            
            
            % Add resistor in parallel to Cap:
            % Add new voltage measurement element to NLnets:
            NLnets(end+1,:) = [strcat(NLnets(switches(j),1),'_R_ON'), NLnets(switches(j),2:end)];
            % Add new voltage measurement (current source) element to NL:
            NL(end+1,:)= [numR,sum(row.*NL(:,2:3)),NL(end,4)+1];
            
            % Set C to be current position:
            % Add new voltage measurement element to NLnets:
            NLnets(switches(j),:) = [strcat(NLnets(switches(j),1),'_C'), NLnets(switches(j),2:end)];
            % Add new voltage measurement (current source) element to NL:
            NL(switches(j),1)= [numC];
            
            % Add voltage source
            row = zeros(size(NL,1),1); % Find element that we are trying to measure the voltage of
            row(end) = 1;
            
            new_node = max(max(NL(:,2:3)))+1; % Set new node number (highest node number + 1)
            node = NL(sum(row.*NL(:,4)),2); % Find node number from element being measured
            
            New_Node = strcat(NLnets(sum(row.*NL(:,4)),2),'_Diode'); % Set new node cell name (highest node number + 1)
            Node = NLnets(sum(row.*NL(:,4)),2); % Find cell name from element being measured
            
            NLnets(sum(row.*NL(:,4)),2) = New_Node; % Adjust element beign measured node (cell net list)
            NLnets(end+1,:) = [strcat(NLnets(switches(j),1),'_VF'),New_Node,Node,NLnets(sum(row.*NL(:,4)),4:end)]; % Create measurement node (cell net list)
            % NLnets(sum(row.*NL(:,4)),1) = strcat(NLnets(switches(j),1),'_R_ON');
            
            NL(sum(row.*NL(:,4)),2)=new_node; % Adjust element beign measured node (numerical net list)
            NL(end+1,:)= [numV,new_node,node,NL(end,4)+1]; % Create measurement node (numerical net list)
            ON_States(ON_index,1) = NLnets(end-1,1);
            ON_index = ON_index+1;
            
            count = count+1;
            
            % BODY DIODE switch is replaced with resistor and forward voltage
            % source have different if statements for diodes and FETs due to
            % FETs diode being opposite in direction given node and sign
            % convention for LTSpice
        end
    end
    
    % if Diode then need to check that it is a diode and that the
    % state is set to zero
    if state(i,j) == 0 && (NL(switches(j),1)==numD)
        row = zeros(size(NL,1),1); % Find element that we are trying to measure the voltage of
        row(switches(j)) = 1;
        
        % Add resistor in parallel to Cap:
        % Add new voltage measurement element to NLnets:
        NLnets(end+1,:) = [strcat(NLnets(switches(j),1),'_R_OFF'), NLnets(switches(j),2:end)];
        % Add new voltage measurement (current source) element to NL:
        NL(end+1,:)= [numR,sum(row.*NL(:,2:3)),NL(end,4)+1];
        
        % Set C to be current position:
        % Add new voltage measurement element to NLnets:
        NLnets(switches(j),:) = [strcat(NLnets(switches(j),1),'_C'), NLnets(switches(j),2:end)];
        % Add new voltage measurement (current source) element to NL:
        NL(switches(j),1)= [numC];
        
        % Add voltage source
        row = zeros(size(NL,1),1); % Find element that we are trying to measure the voltage of
        row(end) = 1;
        
        new_node = max(max(NL(:,2:3)))+1; % Set new node number (highest node number + 1)
        node = NL(sum(row.*NL(:,4)),2); % Find node number from element being measured
        
        New_Node = strcat(NLnets(sum(row.*NL(:,4)),2),'_Diode'); % Set new node cell name (highest node number + 1)
        Node = NLnets(sum(row.*NL(:,4)),2); % Find cell name from element being measured
        
        NLnets(sum(row.*NL(:,4)),2) = New_Node; % Adjust element beign measured node (cell net list)
        NLnets(end+1,:) = [strcat(NLnets(switches(j),1),'_VF'),Node,New_Node,NLnets(sum(row.*NL(:,4)),4:end)]; % Create measurement node (cell net list)
        % NLnets(sum(row.*NL(:,4)),1) = strcat(NLnets(switches(j),1),'_R_OFF');
        
        NL(sum(row.*NL(:,4)),2)=new_node; % Adjust element beign measured node (numerical net list)
        NL(end+1,:)= [numV,node,new_node,NL(end,4)+1]; % Create measurement node (numerical net list)
        ON_States(ON_index,1) = NLnets(end-1,1);
        ON_index = ON_index+1;
        
        
        count = count+1;
        
    end
    
    
    % For power diodes
    
    % Check if state is set to 1 and that it is a diode
    if state(i,j) == 1 && (NL(switches(j),1)==numD)
        
        
        
        row = zeros(size(NL,1),1); % Find element that we are trying to measure the voltage of
        row(switches(j)) = 1;
        
        
        % Add resistor in parallel to Cap:
        % Add new voltage measurement element to NLnets:
        NLnets(end+1,:) = [strcat(NLnets(switches(j),1),'_R_D'), NLnets(switches(j),2:end)];
        % Add new voltage measurement (current source) element to NL:
        NL(end+1,:)= [numR,sum(row.*NL(:,2:3)),NL(end,4)+1];
        
        % Set C to be current position:
        % Add new voltage measurement element to NLnets:
        NLnets(switches(j),:) = [strcat(NLnets(switches(j),1),'_C'), NLnets(switches(j),2:end)];
        % Add new voltage measurement (current source) element to NL:
        NL(switches(j),1)= [numC];
        
        % Add voltage source
        row = zeros(size(NL,1),1); % Find element that we are trying to measure the voltage of
        row(end) = 1;
        
        new_node = max(max(NL(:,2:3)))+1; % Set new node number (highest node number + 1)
        node = NL(sum(row.*NL(:,4)),2); % Find node number from element being measured
        
        New_Node = strcat(NLnets(sum(row.*NL(:,4)),2),'_Diode'); % Set new node cell name (highest node number + 1)
        Node = NLnets(sum(row.*NL(:,4)),2); % Find cell name from element being measured
        
        NLnets(sum(row.*NL(:,4)),2) = New_Node; % Adjust element beign measured node (cell net list)
        NLnets(end+1,:) = [strcat(NLnets(switches(j),1),'_VF'),Node,New_Node,NLnets(sum(row.*NL(:,4)),4:end)]; % Create measurement node (cell net list)
        % NLnets(sum(row.*NL(:,4)),1) = strcat(NLnets(switches(j),1),'_R_D');
        
        NL(sum(row.*NL(:,4)),2)=new_node; % Adjust element beign measured node (numerical net list)
        NL(end+1,:)= [numV,node,new_node,NL(end,4)+1]; % Create measurement node (numerical net list)
        ON_States(ON_index,1) = NLnets(end-1,1);
        ON_index = ON_index+1;
        
        forward_pass(count) = 1;
        count = count+1;
        
        
    end
    
    % For FETS body diode ON
    if j>1 % Without this have an index error because need to check the previous index
        % Have to check if the current state is set to 1 AND if the
        % previous index represents the same state AND if the index
        % corresponds to a FET
        if state(i,j) == 1 && (switches(j-1) == switches(j)) && (NL(switches(j),1)==numM)
            
            
            row = zeros(size(NL,1),1); % Find element that we are trying to measure the voltage of
            row(switches(j)) = 1;
            
            
            % Add resistor in parallel to Cap:
            % Add new voltage measurement element to NLnets:
            NLnets(end+1,:) = [strcat(NLnets(switches(j),1),'_R_D'), NLnets(switches(j),2:end)];
            % Add new voltage measurement (current source) element to NL:
            NL(end+1,:)= [numR,sum(row.*NL(:,2:3)),NL(end,4)+1];
            
            % Set C to be current position:
            % Add new voltage measurement element to NLnets:
            NLnets(switches(j),:) = [strcat(NLnets(switches(j),1),'_C'), NLnets(switches(j),2:end)];
            % Add new voltage measurement (current source) element to NL:
            NL(switches(j),1)= [numC];
            
            % Add voltage source
            row = zeros(size(NL,1),1); % Find element that we are trying to measure the voltage of
            row(end) = 1;
            
            new_node = max(max(NL(:,2:3)))+1; % Set new node number (highest node number + 1)
            node = NL(sum(row.*NL(:,4)),2); % Find node number from element being measured
            
            New_Node = strcat(NLnets(sum(row.*NL(:,4)),2),'_Diode'); % Set new node cell name (highest node number + 1)
            Node = NLnets(sum(row.*NL(:,4)),2); % Find cell name from element being measured
            
            NLnets(sum(row.*NL(:,4)),2) = New_Node; % Adjust element beign measured node (cell net list)
            NLnets(end+1,:) = [strcat(NLnets(switches(j),1),'_VF'),New_Node,Node,NLnets(sum(row.*NL(:,4)),4:end)]; % Create measurement node (cell net list)
            %NLnets(sum(row.*NL(:,4)),1) = strcat(NLnets(switches(j),1),'_R_D');
            
            NL(sum(row.*NL(:,4)),2)=new_node; % Adjust element beign measured node (numerical net list)
            NL(end+1,:)= [numV,new_node,node,NL(end,4)+1]; % Create measurement node (numerical net list)
            ON_States(ON_index,1) = NLnets(end-1,1);
            ON_index = ON_index+1;
            
            forward_pass(count) = 1;
            count = count+1;
            
        end
    end
    
    
    if state(i,j) ~= 1 && state(i,j) ~= 0
        error('Invalid State Detected')
    end
    
end

catch ME
    rethrow(ME);
end


% Update classes
obj.ON_States(:,i)=ON_States;
obj.OFF_States(:,i)=OFF_States;



% Create a variable to multiply B by
% Here this makes an assumption that the order of u is the same as the
% order






end % That's all Folks
