function [] = cutset_loop_num(obj)
%CUTSET_LOOP_NUM Runs after 
%   Detailed explanation goes here

% % % % If this is a simulink file then return without parsing circuit
% % % if isempty(obj.NL)
% % %   return
% % % end
% % % 
% % % %% Find Diodes and Switches
% % % [switches]=obj.findDM;
% % % 
% % % %% Get binary representation of number of states to change R and C for D and M
% % % number_of_states = 2^length(switches);
% % % 
% % % 
% % % state  = 1;
% % % % to find body diodes for DAB
% % % %[state] = obj.bodydiode_correction(switches,state);
% % % 
% % % number_of_states = size(state,1);
% % % 
% % % % Pre-set number of possible Tree and Cotree matrices
% % % %  This is needed to efficiently pass these variables through for each time
% % % %  interval
% % % SortedTree=zeros(2*size(obj.NL,1),5,1);
% % % SortedCoTree=zeros(2*size(obj.NL,1),5,1);
% % % 
% % % % Pre-set number of possible ON and OFF states
% % % %  This is needed to efficiently pass these variables through for each time
% % % %  interval
% % % % These are no longer used in this implementation if there are
% % % % alterations or revision back to a previous version there should be a
% % % % column length of "number_of_states" insead of 1 
% % % obj.ON_States = cell(length(switches),1);
% % % obj.OFF_States = cell(length(switches),1);
% % % 
% % % 
% % % number_of_states = 1;
% % % 
% % % i = 1;
% % % 
% % % [NL,NLnets,forward_pass]=obj.Single_states_D(state,i,switches);

NL = obj.NL;

%% Find incidence matricies
SortedRows = sortrows(NL,1); % sorts rows in prefered tree order
incidence = zeros(max(max(NL(:,2:3))),size(NL,1)); % initialize incidence matrix

% Find Complete Incidence Matrix
for i = 1:1:size(SortedRows,1)
    incidence(SortedRows(i,2),i) = 1;
    incidence(SortedRows(i,3),i) = -1;
end

% Find the Reduced Incidence Matrix
incidence_1 = incidence(2:end,:);
% incidence_2 = incidence(1:end-1,:);

% Reduced Reduced Incidence Matrix to Echelon Form
% REF2 = rref(incidence_2);
REF = rref(incidence_1);

[temp,Ba]=size(incidence_1);
Tree = 0;

%% Find the Normal Tree and CoTree
% Find the first one (1) in each row
% The columns where these ones (1) occur are the braches within the normal
% tree
try
    for i = 1:1:temp
        Tree(i) = find(REF(i,:),1,'first');
    end
catch
    error('Possible Singular Incidence Matrix. Check ground nodes in netlist (No isolation across transformer)')
end


% Find the remaining elements that are not in the normal tree
% The branches are in the Cotree
Ba = 1:Ba;
CoTree = setdiff(Ba,Tree);

%% Orders Tree and CoTree indices correctly: E R G J

numV = 1;
numBV = 2;
numMV = 3;
numC = 4;
numR = 5;
numL = 6;
numMI = 7;
numBI = 8;
numI = 9;


% Adjust Branch Identification numbers from:
% V BV MV  C  R  L MI BI  I
% 1  2  3  4  5  6  7  8  9

% to

% Branch Identification numbers:
% __________Tree__________  _________CoTree__________
% E  E-B  E-M  E-C  E-L  R  G   J-C  J-L J-M  J-B  J
% 1   2    3    4    5   6  7    8    9   10   11  12


% Finds coomponents in tree and reassigns component identifiers:
tosorttree = [SortedRows(Tree(:),:),Tree(:)];
[RowTree,~]=size(tosorttree);
for i = 1:1:RowTree

    switch tosorttree(i,1)
        case numL
            tosorttree(i,1) = 5;
        case numR
            tosorttree(i,1) = 6;
    end

end
SortedTree = sortrows(tosorttree,1);
Tree = SortedTree(:,5);


% Finds coomponents in CoTree and reassigns component identifiers:
tosortcotree = [SortedRows(CoTree(:),:),CoTree(:)];
[RowCoTree,~]=size(tosortcotree);
for i = 1:1:RowCoTree

    switch tosortcotree(i,1)
        case numR
            tosortcotree(i,1) = 7;
        case numC
            tosortcotree(i,1) = 8;
        case numL
            tosortcotree(i,1) = 9;
        case numMI
            tosortcotree(i,1) = 10;
        case numBI
            tosortcotree(i,1) = 11;
        case numI
            tosortcotree(i,1) = 12;
    end

end
SortedCoTree = sortrows(tosortcotree,1);
CoTree = SortedCoTree(:,5);


% For Transformers:
%{
numV = 1;
numBV = 2;
numC = 3;
numR = 4;
numL = 5;
numBI = 6;
numI = 7;

% Branch Identification numbers:
% E  E-B  E-C  E-L  R   G   J-C  J-L  J-B  J
% 1   2    3   4    5   6    7   8     9   10

tosorttree = [SortedRows(Tree(:),:),Tree(:)];
[RowTree,~]=size(tosorttree);
for i = 1:1:RowTree

    switch tosorttree(i,1)
        case numL
            tosorttree(i,1) = 4;
        case numR
            tosorttree(i,1) = 5;
    end

end
SortedTree = sortrows(tosorttree,1);
Tree = SortedTree(:,5);

tosortcotree = [SortedRows(CoTree(:),:),CoTree(:)];
[RowCoTree,~]=size(tosortcotree);
for i = 1:1:RowCoTree

    switch tosortcotree(i,1)
        case numR
            tosortcotree(i,1) = 6;
        case numC
            tosortcotree(i,1) = 7;
        case numL
            tosortcotree(i,1) = 8;
        case numBI
            tosortcotree(i,1) = 9;
        case numI
            tosortcotree(i,1) = 10;
    end

end
SortedCoTree = sortrows(tosortcotree,1);
CoTree = SortedCoTree(:,5);
%}
%{

% Branch Identification numbers:
% E  E-C  E-L  R  G  J-C  J-L  J
% 1   2    3   4  5   6    7   8




tosorttree = [SortedRows(Tree(:),:),Tree(:)];
[RowTree,~]=size(tosorttree);
for i = 1:1:RowTree

    switch tosorttree(i,1)
        case 4
            tosorttree(i,1) = 3;
        case 3
            tosorttree(i,1) = 4;
    end

end
SortedTree = sortrows(tosorttree,1);
Tree = SortedTree(:,5);

tosortcotree = [SortedRows(CoTree(:),:),CoTree(:)];
[RowCoTree,~]=size(tosortcotree);
for i = 1:1:RowCoTree

    switch tosortcotree(i,1)
        case 5
            tosortcotree(i,1) = 8;
        case 3
            tosortcotree(i,1) = 5;
        case 2
            tosortcotree(i,1) = 6;
        case 4
            tosortcotree(i,1) = 7;
    end
end
SortedCoTree = sortrows(tosortcotree,1);
CoTree = SortedCoTree(:,5);

%}

%% Finding A, B, and D Matrix
A_a = [incidence(:,Tree(:)),incidence(:,CoTree(:))]; % Complete incidence matrix
A = [incidence_1(:,Tree(:)),incidence_1(:,CoTree(:))]; % Reduced incidence matrix
A_T = [incidence_1(:,Tree(:))]; % Tree Branches
A_L = [incidence_1(:,Tree(:))]; % Loops
D = inv(A_T)*A; % Fundamental Cutset Matrix
D_L = [REF(:,CoTree(:))];
B_T = -D_L.';
[temp,~] = size(B_T);
B = [B_T,eye(temp)]; % Fundamental Loop Matrix


obj.Cutset = D_L;
obj.SortedTree_cutloop = SortedTree;
obj.SortedCoTree_cutloop = SortedCoTree;


obj.NL = NL;



end

