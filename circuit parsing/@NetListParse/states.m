function [NL,NLnets] = states(obj,state,i,switches)
% STATES finds the different switch positions of the netlist file converter
% and sets the net list for nodeloop by changing open switches to caps and
% on switches to caps in parallel with an on resistance.


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


for j = 1:1:length(switches)
    
    % OPEN switch is replaced with cap
    if j~=length(switches)

        % if a FET then need to check that state and next state is
        % zero AND check that the current and next states correspond
        % to the same switch AND that is is a FET
        
        if (state(i,j) == 0 && state(i,j+1) == 0 && switches(j+1) == switches(j) && (NL(switches(j),1)==numM)) 
            NL(switches(j),1) = numC;
            NLnets(switches(j),:) = [strcat(NLnets(switches(j),1),'_C'), NLnets(switches(j),2:end)];
            OFF_States(OFF_index,1) = NLnets(switches(j),1);
            OFF_index = OFF_index+1;
        end
        
        % CLOSED switch is replaced with cap and resistor
        if (state(i,j) == 1 && j==1) || (state(i,j) == 1 && switches(j+1) == switches(j))
            NL(switches(j),1) = numC;
            [row_NL,~]=size(NL);
            NL(row_NL+1,:) = [numR,NL(switches(j),2:3),row_NL+1];
            NLnets(row_NL+1,:) = [strcat(NLnets(switches(j),1),'_R'), NLnets(switches(j),2:end)];
            NLnets(switches(j),:) = [strcat(NLnets(switches(j),1),'_C'), NLnets(switches(j),2:end)];
            ON_States(ON_index,1) = NLnets(switches(j),1);
            ON_index = ON_index+1;
            
            % BODY DIODE switch is replaced with resistor and forward voltage
            % source have different if statements for diodes and FETs due to
            % FETs diode being opposite in direction given node and sign
            % convention for LTSpice
        end
    end
    
    % if Diode then need to check that it is a diode and that the
    % stae is set to zero
    if state(i,j) == 0 && (NL(switches(j),1)==numD)
        NL(switches(j),1) = numC;
        NLnets(switches(j),:) = [strcat(NLnets(switches(j),1),'_C'), NLnets(switches(j),2:end)];
        OFF_States(OFF_index,1) = NLnets(switches(j),1);
        OFF_index = OFF_index+1;
    end
    
    
    % For power diodes
    
    % Check if state is set to 1 and that it is a diode
    if state(i,j) == 1 && (NL(switches(j),1)==numD)
        NL(switches(j),1) = numR;
        row = zeros(size(NL,1),1); % Find element that we are trying to measure the voltage of
        row(switches(j)) = 1;
        
        new_node = max(max(NL(:,2:3)))+1; % Set new node number (highest node number + 1)
        node = NL(sum(row.*NL(:,4)),2); % Find node number from element being measured
        
        New_Node = strcat(NLnets(sum(row.*NL(:,4)),2),'_Diode'); % Set new node cell name (highest node number + 1)
        Node = NLnets(sum(row.*NL(:,4)),2); % Find cell name from element being measured
        
        NLnets(sum(row.*NL(:,4)),2) = New_Node; % Adjust element beign measured node (cell net list)
        NLnets(end+1,:) = [strcat(NLnets(switches(j),1),'_VF'),Node,New_Node,NLnets(sum(row.*NL(:,4)),4:end)]; % Create measurement node (cell net list)
        
        NL(sum(row.*NL(:,4)),2)=new_node; % Adjust element beign measured node (numerical net list)
        NL(end+1,:)= [numV,node,new_node,NL(end,4)+1]; % Create measurement node (numerical net list)
        ON_States(ON_index,1) = NLnets(switches(j),1);
        ON_index = ON_index+1;
    end
    
    % For FETS body diode
    if j>1 % Without this ahve an index error because need to check the previous index
        % Have to check if the current state is set to 1 AND if the
        % previous index represents the same state AND if the index
        % corresponds to a FET
        if state(i,j) == 1 && (switches(j-1) == switches(j)) && (NL(switches(j),1)==numM)
            
            NL(switches(j),1) = numR;
            row = zeros(size(NL,1),1); % Find element that we are trying to measure the voltage of
            row(switches(j)) = 1;
            
            new_node = max(max(NL(:,2:3)))+1; % Set new node number (highest node number + 1)
            node = NL(sum(row.*NL(:,4)),2); % Find node number from element being measured
            
            New_Node = strcat(NLnets(sum(row.*NL(:,4)),2),'_Diode'); % Set new node cell name (highest node number + 1)
            Node = NLnets(sum(row.*NL(:,4)),2); % Find cell name from element being measured
            
            NLnets(sum(row.*NL(:,4)),2) = New_Node; % Adjust element beign measured node (cell net list)
            NLnets(end+1,:) = [strcat(NLnets(switches(j),1),'_VF'),New_Node,Node,NLnets(sum(row.*NL(:,4)),4:end)]; % Create measurement node (cell net list)
            
            NL(sum(row.*NL(:,4)),2)=new_node; % Adjust element beign measured node (numerical net list)
            NL(end+1,:)= [numV,new_node,node,NL(end,4)+1]; % Create measurement node (numerical net list)
            ON_States(ON_index,1) = NLnets(switches(j),1);
            ON_index = ON_index+1;
            
        end
    end
    
    
    if state(i,j) ~= 1 && state(i,j) ~= 0
        error('Invalid State Detected')
    end
    
end

% Update classes
obj.ON_States(:,i)=ON_States;
obj.OFF_States(:,i)=OFF_States;


end % That's all Folks
