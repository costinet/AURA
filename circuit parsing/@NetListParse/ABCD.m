function [] = ABCD(obj)
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


%% Creat Parse Class

%file = NetListParse();

%% Read in file:
obj.read_file;

%% Find Diodes and Switches

[switches]=obj.findDM;

%% Get binary representation of number of states to change R and C for D and M

number_of_states = 2^length(switches);
bin=de2bi(0:number_of_states-1);
state = bin;


%% Cycle through all possible states
SortedTree=zeros(2*size(obj.NL,1),5,1); 
SortedCoTree=zeros(2*size(obj.NL,1),5,1);
for i = 1:1:1 %number_of_states
    
    [NewNL,NewNLnets]=obj.states(state,i,switches);
    [A(:,:,i),B(:,:,i),C(:,:,i),D(:,:,i),HtempAB(:,:,i),dependsAB(:,:,i),HtempCD(:,:,i),savedCD(:,:,i),StateNamesAB(:,i),StateNamesCD(:,i),OutputNames(:,i),DependentNames(:,i),SortedTree(:,:,i),SortedCoTree(:,:,i)] = obj.nodeloop(NewNL,NewNLnets);
    
    J = 9572839;
end


obj.Asym = A;
obj.Bsym = B;
obj.Csym = C;
obj.Dsym = D;
obj.HtempAB = HtempAB;
obj.HtempCD = HtempCD;
obj.dependsAB = dependsAB;
obj.savedCD = savedCD;
obj.StateNames = StateNamesAB;
obj.OutputNamesCD = StateNamesCD;
obj.DependentNames = DependentNames;
obj.OutputNames = OutputNames;
obj.SortedTree = SortedTree;
obj.SortedCoTree = SortedCoTree;

    J = 5782975892;

end
