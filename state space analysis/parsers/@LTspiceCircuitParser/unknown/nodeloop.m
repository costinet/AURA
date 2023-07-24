function [A,B,C,D,HtempAB,dependsAB,HtempCD,savedCD,StateNamesAB,StateNamesCD,OutputNames,DependentNames,SortedTree1,SortedCoTree1,ConstantNames,OrderedNamesnum] = nodeloop(obj,NL,NLnets)
% nodeloop creates the state matrices A,B,C,D from a specific node input
% matrix
%
% Notes: Needs to be more general, broken up into a few more functions,
%
%


%{
% For the AC fly converter::::
NLnets(11,:) = []
NL(11,:) = []
NL(:,4) = [1:34]
%}


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
    error('ossiPble Singular Incidence Matrix. Check ground nodes in netlist (No isolation across transformer)')
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
A_L = [incidence_1(:,CoTree(:))]; % Link Branches
D = inv(A_T)*A; % Fundamental Cutset Matrix
D_L = [REF(:,CoTree(:))]; % Links Matrix
B_T = -D_L.';
[temp,~] = size(B_T);
B = [B_T,eye(temp)]; % Fundamental Loop Matrix

[row,col]=size(D_L);

%% Finding size of Z and G Matrices

% New Index Values
% Branch Identification numbers:
% __________Tree__________  _________CoTree__________
% E  E-B  E-M  E-C  E-L  R  G   J-C  J-L J-M  J-B  J
% 1   2    3    4    5   6  7    8    9   10   11  12

numE = 1;
numEB = 2;
numEM = 3;
numEC = 4;
numEL = 5;
numR = 6;
numG = 7;
numJC = 8;
numJL = 9;
numJM = 10;
numJB = 11;
numJ = 12;


hey = find(SortedTree(:,1)==numR,1); % find first index of resistor (R) in tree
if isempty(hey) % if there are no R in tree then set to zero
    R = row+1;
    Z_R = 0;
    Y_R = 0;
else % if there are resistors in tree
    % Create symbolic values for all resistor names:
    R = hey;
    syms(NLnets(SortedTree(hey:end,4),1));
    syms 'Z_R';
    % Create symbolic resistor matrix for Tree:
    for i = 1:1:row-R+1
        Z_R(i,i) = NLnets(SortedTree(R+i-1,4),1);
    end
    % Find inverse of symbolic resistor matrix for Tree:
    Y_R = inv(Z_R);
end

ya = find(SortedCoTree(:,1)==numG,1,'last'); % find last index of resistor (R) in cotree
if isempty(ya) % if there are no resistors in cotree then set to zero
    G = 0;
    Y_G = 0;
    Z_G = 0;
else % if there are resistors in cotree
    % Create symbolic values for all resistor names:
    G = ya;
    syms(NLnets(SortedCoTree(1:ya,4),1));
    syms 'Z_G';
    % Create symbolic resistor matrix for CoTree:
    for i = 1:1:G
        Z_G(i,i) = NLnets(SortedCoTree(i,4),1);
    end
    % Find inverse of symbolic resistor matrix for CoTree:
    Y_G = inv(Z_G);
end

%% Parsing D_L Matrix
D_EG = D_L(1:R-1,1:G);
D_EJ = D_L(1:R-1,G+1:col);
D_RG = D_L(R:row,1:G);
D_RJ = D_L(R:row,G+1:col);

% D is the fundamental cutset matrix
%
%       E     R      G     J
%      _                     _
%     |1_EE  0_ER | D_EG  D_EJ|
% D = |           |           |
%     |0_RE  1_RR | D_RG  D_RJ|
%      -                     -

% Initalize variables
Flag_EG = 0;
Flag_EJ = 0;
Flag_RG = 0;
Flag_RJ = 0;

% Check to see if any D matrices are empty and replace with zero and set
% flag
if isempty(D_EG)
    D_EG = 0;
    Flag_EG = 1;
end
if isempty(D_EJ)
    D_EJ = 0;
    Flag_EJ = 1;
end
if isempty(D_RG)
    D_RG = 0;
    Flag_RG = 1;
end
if isempty(D_RJ)
    D_RJ = 0;
    Flag_RJ = 1;
end

%% Solve for Z and Y Matrices
Z = Z_G+(D_RG.')*Z_R*D_RG;
Y = Y_R+D_RG*Y_G*(D_RG.');

% Solve for the inverse of Z and Y Matrices
if Z == 0
    inv_Z = 0;
else
    inv_Z = inv(Z);
end

if Y == 0
    inv_Y = 0;
else
    inv_Y = inv(Y);
end

%% Find H Sub Matrix

% Add more cases as needed for now


if Flag_RJ == 0 && Flag_RG == 0 && Flag_EG == 0 && Flag_EJ == 0
    % Normal Case
    H_EE = -(D_EG*inv_Z*(D_EG.'));
    H_EJ = D_EG*inv_Z*(D_RG.')*Z_R*D_RJ-D_EJ;
    H_JE = (D_EJ.')-((D_RJ.')*inv_Y*D_RG*Y_G*(D_EG.'));
    H_JJ = -((D_RJ.')*inv_Y*D_RJ);
    almost_H = [H_EE,H_EJ;H_JE,H_JJ];
else
    if Flag_RJ == 1 && Flag_RG == 1 && Flag_EG == 0 && Flag_EJ == 0
        % Correct H_JJ if there are no resistors in the tree
        H_EE = -(D_EG*inv_Z*(D_EG.'));
        H_EJ = -D_EJ;
        H_JE = (D_EJ.');
        [~,H_EJ_col]=size(H_EJ);
        [H_JE_row,~]=size(H_JE);
        H_JJ = zeros(H_JE_row,H_EJ_col);
        almost_H = [H_EE,H_EJ;H_JE,H_JJ];
    else
        if  Flag_EG == 1 && Flag_RG == 1 && Flag_RJ == 0 && Flag_EJ == 0
            % Correct H_EE if there are no resistors in the cotree
            H_JJ = -((D_RJ.')*inv_Y*D_RJ);
            H_EJ = -D_EJ;
            H_JE = (D_EJ.');
            [H_EJ_row,~]=size(H_EJ);
            [~,H_JE_col]=size(H_JE);
            H_EE = zeros(H_EJ_row,H_JE_col);
            almost_H = [H_EE,H_EJ;H_JE,H_JJ];
        else
            if Flag_RJ == 1 && Flag_RG == 1 && Flag_EG == 1 && Flag_EJ == 0
                H_EJ = -D_EJ;
                H_JE = (D_EJ.');
                [H_EJ_row,H_EJ_col]=size(H_EJ);
                [H_JE_row,H_JE_col]=size(H_JE);
                H_JJ = zeros(H_JE_row,H_EJ_col);
                H_EE = zeros(H_EJ_row,H_JE_col);
                almost_H = [H_EE,H_EJ;H_JE,H_JJ];
            else
                error('Unknown combination of elements')
            end
        end
    end
end

%% Find H and s Matrix

% H matrix including voltage sources
almost_H = [H_EE,H_EJ;H_JE,H_JJ];
[H,s]=obj.hybridparse(almost_H,SortedTree,SortedCoTree);

% Eventual Function to find outputs
[A,B,C,D,HtempAB,dependsAB,StateNamesAB,OutputNames,DependentNames,ConstantNames,OrderedNamesnum]=obj.loopfixAB(H,s,NLnets,SortedTree,SortedCoTree);

[C,D,HtempCD,savedCD,StateNamesCD]=obj.loopfixCD(A,B,C,D,H,s,NLnets,SortedTree,SortedCoTree);

SortedTree1=zeros(2*size(obj.NL,1),5,1);
SortedCoTree1=zeros(2*size(obj.NL,1),5,1);
SortedTree1(1:size(SortedTree,1),1:5) = SortedTree;
SortedCoTree1(1:size(SortedCoTree,1),1:5) = SortedCoTree;
end

%{
% find and remove voltage sources to find h and s matrix
lastvsource = find(SortedTree(:,1)==1,1,'last');
firstisource = find(SortedCoTree(:,1)==8,1,'first');

if isempty(lastvsource)
    H = almost_H;
    s = 0;
else
    H = almost_H(lastvsource+1:end,lastvsource+1:end);
    s = almost_H(lastvsource+1:end,1:lastvsource);
end
%}
%{
if ~isempty(firstisource)
    H = H(1:firstisource-1,1:firstisource-1

end
%}
%{
% Create identity matrix along with H and almost H matrix in order to find
% correct orientation of state and output matrix
[H_row,H_col] = size(H); % Get size of Hybrid Matrix without s
[almostH_row,~] = size(almost_H); % Get size of full Hybrid Matrix
Htemp = [eye(H_row),H,s]; % Create Mx = Ax + Bu form
almostHtemp = [eye(almostH_row),almost_H]; % Create Mx = Ax + Bu form

temps = [SortedTree(:,:);SortedCoTree(:,:)];
cir = temps(:,1)~=numR & temps(:,1)~=numG; % find position of elements in almost H
cir_state = find(temps(:,1)~=numR & temps(:,1)~=numG & temps(:,1)~=numE & temps(:,1)~=numJ & temps(:,1)~=numEB & temps(:,1)~=numJB); % find postitin of elements in H

OrderedNamesnum = temps(cir_state,4); % Find the index for the output state names
OutputNames = NLnets(OrderedNamesnum,1); % Find the list of output names
DependentNames = {};

OrderedNameselement = temps(cir_state,1);
loop = length(OrderedNameselement)+1; % Set while loop index

% break somewhere in here to separate different solves for AB matrix and CD
% matrix

% For loop to switch i and v in the hybrid matrix to ensure all caps have
% an output current and all inductors have an output voltage:
for i = 1:1:length(OrderedNameselement)
 if OrderedNameselement(i)==numEL || OrderedNameselement(i)==numJC
        [H_row2,~]=size(Htemp);
        move=Htemp(:,i);
        move2 = Htemp(:,i+H_row2);
        Htemp(:,i)=-move2;
        Htemp(:,i+H_row2) = -move;
 end
end

% Pluggs in all the C and L values to from M*x_dot = Ax+Bu equations
% Basically if the number exits in the M matrix then the entire row gets
% devided by that elements either L or C
% The position for that element however does not since it forms either
% di/dt or dv/dt

for i = 1:1:length(OrderedNameselement)
    for j = 1:1:length(OrderedNameselement)
        if Htemp(i,j)~=0
            Htemp(i,:)=Htemp(i,:)./OutputNames(j);
            Htemp(i,j) = Htemp(i,j).* OutputNames(j);
        end
    end
end

i = 1;
j = 0;
k = 0;

while i<loop
%depends = [];
%for i = 1:1:length(OrderedNameselement)
    % need to switch if caps were in cotree or inductors in tree
    k = k+1;
    if OrderedNameselement(k)==numEL || OrderedNameselement(k)==numJC

        DependentNames(j+1,:) = OutputNames(i); % Assigns name to dependent row status
        OutputNames(i) = []; % Deletes name from Output names list

        j = j+1; % J is number of dependent elements

        % Set up so each depends row is a linear combination of the
        % independent state remaining:
        depends(j,:)=Htemp(i,H_row2+1:2*H_row2);
        depends(j,i) = 0;

        Htemp(1:H_row2,1:H_row2) = Htemp(1:H_row2,1:H_row2) + repmat(depends(j,:),H_row2,1).*Htemp(:,i);
        Htemp(i,:)=[];
        Htemp(:,i+H_row2) = [];
        Htemp(:,i)=[];
        depends(:,i)=[];

        i = i-1;
        loop = loop - 1;
        J = 5829758239;
        [H_row2,~] = size(Htemp);

    end
    i= i+1;
% end
end

% Need to format output how i want it



[H_row2,~] = size(Htemp);
H_row2 = H_row2+1;

% Multiply columns times the state variables
syms(OutputNames);

% for i = 1:1:length(OutputNames)



%Htemp(:,1:H_row2-1) = OutputNames(:)'.*Htemp(:,1:H_row2-1);

% end

% Cannot do this for large matrix (will have to do after eval)
Htemp = rref(Htemp);

OutputHtemp = Htemp;

% % statedepends = [];
% % statedependsconst = [];
% % outstatedepends = [];
% % outstatedependsconst = [];

for i = 1:1:j
    dependstate = depends(i,:)'.*Htemp(:,H_row2:2*(H_row2-1));
    dependsconst = depends(i,:)'.*Htemp(:,(2*(H_row2-1))+1:end);
    statedependsconst(i,:) = sum (dependsconst);
    statedepends(i,:)=sum(dependstate);
    % add all columns of matrix to get equation that goes in state equation
end

A = [Htemp(:,H_row2:2*(H_row2-1)),zeros(H_row2-1,j);statedepends(:,:),zeros(j,j)];

% fix B
B = [Htemp(:,(2*(H_row2-1))+1:end);statedependsconst];


OutputHtemp = rref(OutputHtemp);

StateNames = [OutputNames;DependentNames];


for i = 1:1:j
    outdependstate = depends(i,:)'.*Htemp(:,H_row2:2*(H_row2-1));
    outdependsconst = depends(i,:)'.*Htemp(:,(2*(H_row2-1))+1:end);
    outstatedependsconst(i,:) = sum(outdependsconst).*DependentNames(i);
    outstatedepends(i,:)=sum(outdependstate).*DependentNames(i);
    % add all columns of matrix to get equation that goes in state equation
end

C = [OutputHtemp(:,H_row2:2*(H_row2-1)).*OutputNames,zeros(H_row2-1,j);outstatedepends(:,:),zeros(j,j)];

D = [OutputHtemp(:,(2*(H_row2-1))+1:end).*OutputNames;outstatedependsconst];

%}
