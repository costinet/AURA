function [A,B,C,D,HtempAB,dependsAB,HtempCD,savedCD,StateNamesAB,StateNamesCD,OutputNames,DependentNames,SortedTree1,SortedCoTree1,ConstantNames,OrderedNamesnum] = nodeloop_num(obj,NL,NLnets)
% nodeloop creates the state matrices A,B,C,D from a specific node input
% matrix
%
% Notes: Needs to be more general, broken up into a few more functions,
%
%

%% Sort numeric index values 
% This section matches the user given values to all variables from the
% topology level to allow numerical solving through a lookup table

% longcellarray = obj.Component_Values(:,1);
% shortcellarray  = NLnets(:,1);
% for ind =1:length(shortcellarray)
% IndexC = strfind(longcellarray, cell2mat(shortcellarray(ind)));
% Index{1,ind} = find(not(cellfun('isempty', IndexC)));
% end


shortcellarray = obj.Component_Values(:,1);
longcellarray  = NLnets(:,1);
index  = zeros(length(longcellarray),1);
for ind =1:length(shortcellarray)
    for ind2 = 1:length(longcellarray)
        truth = strcmp(longcellarray{ind2}, shortcellarray{ind});
        if truth
            index(ind2) = ind;
            break
        end
    end
end
obj.index = index;



 D_L = obj.Cutset;
 SortedTree = obj.SortedTree_cutloop;
 SortedCoTree = obj.SortedCoTree_cutloop;

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
    %syms(NLnets(SortedTree(hey:end,4),1));
    %syms 'Z_R';
    % Create numeric resistor matrix for Tree:
    for i = 1:1:row-R+1
        Z_R(i,i) = obj.Component_Values{index(SortedTree(R+i-1,4)),2};
        % Z_R(i,i) = NLnets(SortedTree(R+i-1,4),1); % What the
        % symbolic vesion was
    end
    % Find inverse of symbolic resistor matrix for Tree:
    Y_R = inv(Z_R);
end

ya = find(SortedCoTree(:,1)==numG,1,'last'); % find last index of resistor (R) in tree
if isempty(ya) % if there are no resistors in cotree then set to zero
    G = 0;
    Y_G = 0;
    Z_G = 0;
else % if there are resistors in cotree
    % Create symbolic values for all resistor names:
    G = ya;
    % syms(NLnets(SortedCoTree(1:ya,4),1));
    % syms 'Z_G';
    % Create symbolic resistor matrix for CoTree:
    for i = 1:1:G
        Z_G(i,i) = obj.Component_Values{index(SortedCoTree(i,4)),2};
        % Z_G(i,i) = NLnets(SortedCoTree(i,4),1);
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
[A,B,C,D,HtempAB,dependsAB,StateNamesAB,OutputNames,DependentNames,ConstantNames,OrderedNamesnum]=obj.loopfixAB_num(H,s,NLnets,SortedTree,SortedCoTree);

[C,D,HtempCD,savedCD,StateNamesCD]=obj.loopfixCD_num(A,B,C,D,H,s,NLnets,SortedTree,SortedCoTree);

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
