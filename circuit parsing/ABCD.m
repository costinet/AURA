function [A,B,C,D,StateNames] = ABCD(filename)
% ABCD creates takes a NETlist file from LTSpice and creates the
% associated ABCD matrices.

%{
% fileID = -1;
% errmsg = '';
% while fileID < 0 
%    disp(errmsg);
%    filename = input('Open file: ', 's');
%    [fileID,errmsg] = fopen(filename);in
% end
%}

%% Read in file:
[NLwhole,NLnets,NL,K]=read_file(filename);


%% Find Diodes


numV = 1;
numBV = 2;
numC = 3;
numR = 4;
numL = 5;
numBI = 6;
numI = 7;
numD = 8;
numM = 9;


switches = [];
Diodes = [];
for i = 1:1:length(NL)
    if NL(i,1) == numD 
        switches(end+1) = i;
        Diodes(end+1,:)  = NL(i,2:3);
    end
    if NL(i,1) == numM
        switches(end+1) = i;
        Diodes(end+1,:)  = [NL(i,3),NL(i,2)];
    end
end

%{ 
'Diodes' contains the position for all of the diodes in the circuit.
These are sorted by:
 Anodes in the 1st column
 Cathodes in the 2nd column
Therefore for every row, there can not be a positive voltage from the node
in the 1st column to the node in the second column
%}

%% Get binary representation of number of states to change R and C for D and M

number_of_states = 2^length(switches);
bin=de2bi(0:number_of_states-1);
state = bin;

% state = bin+2;
% ST = [];

% for i = 1:1:number_of_states
%     for j = 1:1:length(switches)
%         NL(switches(j),1) = state(i,j);
%         % each row is a different state
%     end
%      nodeloop(NL,NLnets);
%     %[NL1(:,:,i),ST(:,:,i)] = Tree(NL,NLnets);
% 
% end

%% Cycle through all possible states

for i = 1:1:number_of_states

    [NewNL,NewNLnets]=states(NL,NLnets,state,i,switches);
    %circuitplot(NewNL,NewNLnets);
    [A(:,:,i),B(:,:,i),C(:,:,i),D(:,:,i),StateNames] = nodeloop(NewNL,NewNLnets,K);

    J = 9572839;
end
    