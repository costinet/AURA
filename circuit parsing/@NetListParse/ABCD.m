function [] = ABCD(obj)
% ABCD creates takes a NETlist file from LTSpice and creates the associated
% ABCD matrices or the associated matrices need to get ABCD if a symbolic
% formulation cannot be solved in a descent amount of time.
%
% A NetListParse class is required with the associated filename defined
%
% filename contains a string of the netlist filename from LTSpice Example:
% boost.net


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
obj.read_file;

% If this is a simulink file then return without parsing circuit
if isempty(obj.NL)
  return
end

%% Find Diodes and Switches
[switches]=obj.findDM;

%% Get binary representation of number of states to change R and C for D and M
number_of_states = 2^length(switches);
bin=de2bi(0:number_of_states-1);
state = bin;

% Pre-set number of possible Tree and Cotree matrices
%  This is needed to efficiently pass these variables through for each time
%  interval
SortedTree=zeros(2*size(obj.NL,1),5,1);
SortedCoTree=zeros(2*size(obj.NL,1),5,1);

% Pre-set number of possible ON and OFF states
%  This is needed to efficiently pass these variables through for each time
%  interval
obj.ON_States = cell(length(switches),number_of_states);
obj.OFF_States = cell(length(switches),number_of_states);

%% Cycle through all possible states

for i = 1:1:number_of_states

    [NewNL,NewNLnets]=obj.states(state,i,switches);
    [A(:,:,i),B(:,:,i),C(:,:,i),D(:,:,i),HtempAB(:,:,i),dependsAB(:,:,i),HtempCD(:,:,i),savedCD(:,:,i),StateNamesAB(:,i),StateNamesCD(:,i),OutputNames(:,i),DependentNames(:,i),SortedTree(:,:,i),SortedCoTree(:,:,i)] = obj.nodeloop(NewNL,NewNLnets);

    J = 9572839;
end

% Update class
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

end % That's all Folks
