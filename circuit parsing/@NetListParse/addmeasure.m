function [] = addmeasure(obj)
%addmeasure(NL,NLnets) adds measurement nodes to NL and NLnets
%   addmeasure takes the numerical net list (NL) and the cell netlist
%   (NLnets) and produces a new NL and NLnets based on the desired voltage
%   and current measurements for the circuit. Currently the desired voltage
%   and current measurements are provided in the function by changing the
%   Voltage and Current variables. Each variable should be a cell array
%   containing the cell netlist identifier for the nodes wanted to 


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

%% Set index Values

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

%% Set Voltage and Current Nodes to add

% Change Voltage and Current based on desired output measurements (C and D
% matricies). Voltage and Current should be of type Cell 
% Example:
% Voltage = {'V1'
%     'M1'
%     'L3'};
% 
% Current = {'C1'
%     'D2'
%     'R3'};


NL = obj.NL;
NLnets = obj.NLnets;
Voltage = obj.Meas_Voltage;
Current = obj.Meas_Current;

%% Set Voltage Measurement Elements

% Voltage Measurement elements are currents sources where i=0 placed in
% parallel to the desired element to measure. This allows the nodes to
% have no impact on the circuit or the formulation of trees and cotrees.

for i = 1:1:length(Voltage)
    
    % Find element that we are trying to measure voltage of:
    row = strcmp(NLnets(:,1),Voltage(i));
    
    if sum(row) == 1
        
        % Add new voltage measurement element to NLnets:
        NLnets(end+1,:) = [strcat(Voltage(i),' V'), NLnets(sum(row.*NL(:,4)),2:end)];
        
        % Add new voltage measurement (current source) element to NL:
        NL(end+1,:)= [numMI,sum(row.*NL(:,2:3)),NL(end,4)+1];
        
    else
        if sum(row) == 0
            warning('Element %s does not exist in the netlist provided \nVoltage measurement %s not implemented into net list',Voltage{i},Voltage{i})
        end
        if sum(row) > 1
            warning('Multiple element %ss exist in the netlist provided \nVoltage measurement %s not implemented into net list',Voltage{i},Voltage{i})
        end
    end
end
%% Set Current Measurement Elements

% Current Measurement elements are voltage sources where v=0 placed in
% series to the desired element to measure. This allows the nodes to
% have no impact on the circuit or the formulation of trees and cotrees.

% To place an additional element in series a new node must be added to the
% circuit. 

for i = 1:1:length(Current)
    row = strcmp(NLnets(:,1),Current(i)); % Find element that we are trying to measure the voltage of
    
    if sum(row) == 1
        
        new_node = max(max(NL(:,2:3)))+1; % Set new node number (highest node number + 1)
        node = NL(sum(row.*NL(:,4)),2); % Find node number from element being measured
        
        New_Node = strcat(NLnets(sum(row.*NL(:,4)),2),'_Meas'); % Set new node cell name (highest node number + 1)
        Node = NLnets(sum(row.*NL(:,4)),2); % Find cell name from element being measured
        
        NLnets(sum(row.*NL(:,4)),2) = New_Node; % Adjust element beign measured node (cell net list)
        NLnets(end+1,:) = [strcat(Current(i),' A'),Node,New_Node,NLnets(sum(row.*NL(:,4)),4:end)]; % Create measurement node (cell net list)
        
        NL(sum(row.*NL(:,4)),2)=new_node; % Adjust element beign measured node (numerical net list)
        NL(end+1,:)= [numMV,node,new_node,NL(end,4)+1]; % Create measurement node (numerical net list)
        
    else
        if sum(row) == 0
            warning('Element %s does not exist in the netlist provided \nCurrent measurement %s not implemented into net list',Current{i},Current{i})
        end
        if sum(row) > 1
            warning('Multiple element %ss exist in the netlist provided \nCurrent measurement %s not implemented into net list',Current{i},Current{i})
        end
    end

end

obj.NL = NL;
obj.NLnets = NLnets;

end % That's all Folks
