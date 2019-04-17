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

[state] = obj.bodydiode_correction(switches,state);

number_of_states = size(state,1);

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

%[THE_LIST]=obj.states_find(state,switches);

ON = [1 0];
OFF = [0 0];

Binary_for_DAB = [
    OFF ON ON OFF OFF ON ON OFF % Reverse power
    OFF OFF OFF OFF OFF ON ON OFF % primary sw
    ON OFF OFF ON OFF ON ON OFF % phase shift
    ON OFF OFF ON OFF OFF OFF OFF % secondary sw
    ON OFF OFF ON ON OFF OFF ON % POWER
    OFF OFF OFF OFF ON OFF OFF ON %primary sw
    OFF ON ON OFF ON OFF OFF ON % phase shift
    OFF ON ON OFF OFF OFF OFF OFF]; % secondary sw

The_Codex = []; % Initialize matrix
    for i = 1:1:size(Binary_for_DAB,1)
        Binary_for_DAB_test = Binary_for_DAB(i,:);
        The_Codex(i) = find(sum(repmat(Binary_for_DAB_test,size(state,1),1)==state,2)==size(state,2)==1);
       
    end


Combinations = The_Codex;

% From THE_LIST, indicate by indicies which state you are wanting to
% model. For eample, if you do not want o
%Combinations = [2 1 4 1];  
[b,m1,~] = unique(Combinations,'first');
[~,d1] =sort(m1);
Combinations = b(d1);


[Combinations] = obj.BDCheck(Combinations,state,switches);

state_old = state;
state = state(Combinations,:);

number_of_states = size(Combinations,2);

%% Cycle through all possible states

for i = 1:1:number_of_states

    [NewNL,NewNLnets,forward_pass]=obj.states(state,i,switches);
    [A(:,:,i),B(:,:,i),C(:,:,i),D(:,:,i),HtempAB(:,:,i),dependsAB(:,:,i),HtempCD(:,:,i),savedCD(:,:,i),StateNamesAB(:,i),StateNamesCD(:,i),OutputNames(:,i),DependentNames(:,i),SortedTree(:,:,i),SortedCoTree(:,:,i),ConstantNames(:,i),OrderedNamesnum(:,i)] = obj.nodeloop(NewNL,NewNLnets);
    
    key = [ones(size(B(:,:,i),2)-size(forward_pass,1)-sum(SortedCoTree(:,1)==12),1);forward_pass;ones(sum(SortedCoTree(:,1)==12),1)]'; % Fit the forward voltage for body diodes in between other voltage sources and current sources to achive corret orientation in B
    B(:,:,i)=B(:,:,i).*repmat(key,size(B(:,:,i),1),1); % Multiply the key times B to get rid of the columns of B where there is not a diode on during ith state
    
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
obj.Codex = The_Codex;

end % That's all Folks
