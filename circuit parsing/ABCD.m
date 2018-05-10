function [A,B,C,D,NLnets,StateNamesAB,StateNamesCD] = ABCD(filename)
% ABCD creates takes a NETlist file from LTSpice and creates the associated
% ABCD matrices.
% 
% 
% filename contians a string of the netlist filename from LTSpice Example:
% boost.net
% 
% A,B,C,D are the coefficents for the state space equation:
%     x(n+1) = A*x(n)+B*u y = C*x(n)+D*u
% 
% NLnets contains a cell array of strings desribing all of the elements in
% the circuit
%
% StateName is the name of the output variables y from the state equations
% 
% 'Diodes' contains the position for all of the diodes in the circuit.
% These are sorted by:
%  Anodes in the 1st column Cathodes in the 2nd column
% Therefore for every row, there can not be a positive voltage from the
% node in the 1st column to the node in the second column
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


%% Read in file:
[NLwhole,NLnets,NL,K]=read_file(filename);


%% Set 

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

switches = [];
Diodes = [];
for i = 1:1:size(NL,1)
    if NL(i,1) == numD 
        switches(end+1) = i;
        Diodes(end+1,:)  = NL(i,2:3);
    end
    if NL(i,1) == numM
        switches(end+1) = i;
        Diodes(end+1,:)  = [NL(i,3),NL(i,2)];
    end
end

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
    [A(:,:,i),B(:,:,i),C(:,:,i),D(:,:,i),StateNamesAB(:,i),StateNamesCD(:,i)] = nodeloop(NewNL,NewNLnets,K);

    J = 9572839;
end

    J = 5782975892;
    
end