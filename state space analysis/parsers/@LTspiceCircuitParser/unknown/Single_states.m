function [NL,NLnets,forward_pass] = Single_states(obj,state,i,switches)
% Single_states finds the different switch positions of the netlist file converter
% and sets the net list for nodeloop by changing open switches to caps and
% on switches to caps in parallel with an on resistance.

% This is different from states because it will add '_R' to resistors
% instead of distinguising ON, OFF, and Diode States 

% This code will result in a single netlist to be parsed symboliclly

%     _   _   _  ____    _
%    / \ | | | |/ _  |  / \
%   / _ \| | | | (_| | / _ \
%  / ___ | |_| |> _  |/ ___ \
% /_/   \_\___//_/ |_/_/   \_\

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

for j = 1:1:length(switches)
    
    % OPEN switch is replaced with cap
    % if a FET then need to check that state and next state is
    % zero AND check that the current and next states correspond
    % to the same switch AND that is is a FET
    
    % if Diode then need to check that it is a diode and that the
    % state is set to zero
    if  (NL(switches(j),1)==numD)
        row = zeros(size(NL,1),1); % Find element that we are trying to measure the voltage of
        row(switches(j)) = 1;
        
        % Add resistor in parallel to Cap:
        % Add new voltage measurement element to NLnets:
        NLnets(end+1,:) = [strcat(NLnets(switches(j),1),'_R'), NLnets(switches(j),2:end)];
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
    
    if j==length(switches)
        
        break
    end
    if  (switches(j+1) == switches(j)) && (NL(switches(j),1)==numM)
        row = zeros(size(NL,1),1); % Find element that we are trying to measure the voltage of
        row(switches(j)) = 1;
        
        
        % Add resistor in parallel to Cap:
        % Add new voltage measurement element to NLnets:
        NLnets(end+1,:) = [strcat(NLnets(switches(j),1),'_R'), NLnets(switches(j),2:end)];
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
        
        %%%%%%%%%%%%%%% New code to put a resistor in series with the FETs
        %{
        
        if switches(j) == 11
        
        new_resist_node = max(max(NL(:,2:3)))+1; % Set new node number (highest node number + 1)
        node_resist = NL(switches(j),2); % Find node number from element being measured
        
        New_Resistor_Node = strcat(NLnets(switches(j),2),'_Resist'); % Set new node cell name (highest node number + 1)
        Node_Resist = NLnets(switches(j),2); % Find cell name from element being measured
        
        NLnets(switches(j),2) = New_Resistor_Node; % Adjust element beign measured node (cell net list)
        NLnets(end+1,:) = [strcat(NLnets(switches(j),1),'_Resist'),Node_Resist,New_Resistor_Node,NLnets(switches(j),4:end)]; % Create measurement node (cell net list)
        
        NL(switches(j),2)=new_resist_node; % Adjust element beign measured node (numerical net list)
        NL(end+1,:)= [numR,node_resist,new_resist_node,NL(end,4)+1]; % Create measurement node (numerical net list)
        
        end
        %}
        %%%%%%%%%%%%%%%
        
        
        ON_States(ON_index,1) = NLnets(end-1,1);
        ON_index = ON_index+1;
        count = count+1;
        j = j+1;
    end
    
    
   % Eliminated state because it is too big when looking at large circuits a lot of circuits 
   % if state(i,j) ~= 1 && state(i,j) ~= 0
   %     error('Invalid State Detected')
   % end
    
end


% Update classes - These don't matter in this implementation
obj.ON_States(:,i)=ON_States;
obj.OFF_States(:,i)=OFF_States;


end % That's all Folks
