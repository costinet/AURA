function [switches] = findDM(obj)
% findDM find the diodes and switches in a circuit 
% findDM should be run after read_file.m 

% Diodes contains the position for all of the diodes in the circuit.
% These are sorted by:
%  Anodes in the 1st column Cathodes in the 2nd column
% Therefore for every row, there can not be a positive voltage from the
% node in the 1st column to the node in the second column

% Switches contains the row location of all switching elements in the
% circuit using obj.NL as a reference

%     _   _   _  ____    _    
%    / \ | | | |/ _  |  / \   
%   / _ \| | | | (_| | / _ \  
%  / ___ | |_| |> _  |/ ___ \ 
% /_/   \_\___//_/ |_/_/   \_\

%% Set Variables

% Element indexes
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

%% Find Diodes and Switches

% Get net lists from class
NL = obj.NL;
NLnets = obj.NLnets;

% Initialize variables
switches = [];
Diodes = [];

% Iterate through netlist to find diodes and switches
for i = 1:1:size(NL,1)
    if NL(i,1) == numD
        switches(end+1) = i;
        Diodes(end+1,:)  = NL(i,2:3);
    end
    if NL(i,1) == numM
        switches(end+1) = i;
        Diodes(end+1,:)  = [NL(i,3),NL(i,2)];
        switches(end+1) = i;
    end
end

% Update class (I'm don't think these hold any value after this)
obj.Diodes = Diodes;
obj.Switches = switches;

end % That's all Folks
