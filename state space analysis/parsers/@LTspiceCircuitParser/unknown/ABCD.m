function [A,B,C,D] = ABCD(obj)
% ABCD creates takes a NETlist file from LTSpice and creates the associated
% ABCD matrices or the associated matrices need to get ABCD if a symbolic
% formulation cannot be solved in a descent amount of time.
%
% A NetListParse class is required with the associated filename defined
%
% filename contains a string of the netlist filename from LTSpice Example:
% boost.net


%     _   _   _  ____    _    
%    / \ | | | |/ _  |  / \   
%   / _ \| | | | (_| | / _ \  
%  / ___ | |_| |> _  |/ ___ \ 
% /_/   \_\___//_/ |_/_/   \_\


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


state  = 1;
% to find body diodes for DAB
%[state] = obj.bodydiode_correction(switches,state);

number_of_states = size(state,1);

% Pre-set number of possible Tree and Cotree matrices
%  This is needed to efficiently pass these variables through for each time
%  interval
SortedTree=zeros(2*size(obj.NL,1),5,1);
SortedCoTree=zeros(2*size(obj.NL,1),5,1);

% Pre-set number of possible ON and OFF states
%  This is needed to efficiently pass these variables through for each time
%  interval
% These are no longer used in this implementation if there are
% alterations or revision back to a previous version there should be a
% column length of "number_of_states" insead of 1 
obj.ON_States = cell(length(switches),1);
obj.OFF_States = cell(length(switches),1);


number_of_states = 1;


%% Cycle through all possible states

for i = 1:1:number_of_states

    [NewNL,NewNLnets,forward_pass]=obj.Single_states(state,i,switches);
   
    % symbolic
    [A(:,:,i),B(:,:,i),C(:,:,i),D(:,:,i),HtempAB(:,:,i),dependsAB(:,:,i),HtempCD(:,:,i),savedCD(:,:,i),StateNamesAB(:,i),StateNamesCD(:,i),OutputNames(:,i),DependentNames(:,i),SortedTree(:,:,i),SortedCoTree(:,:,i),ConstantNames(:,i),OrderedNamesnum(:,i)] = obj.nodeloop(NewNL,NewNLnets);
    
    % numeric
    %[A(:,:,i),B(:,:,i),C(:,:,i),D(:,:,i),HtempAB(:,:,i),dependsAB(:,:,i),HtempCD(:,:,i),savedCD(:,:,i),StateNamesAB(:,i),StateNamesCD(:,i),OutputNames(:,i),DependentNames(:,i),SortedTree(:,:,i),SortedCoTree(:,:,i),ConstantNames(:,i),OrderedNamesnum(:,i)] = obj.nodeloop_num(NewNL,NewNLnets);
   
    J = 98631;
end

% openvar('variable_name')

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
obj.ConstantNames = ConstantNames;
obj.OrderedNamesnum = OrderedNamesnum;
%obj.Codex = The_Codex;

end % That's all Folks
