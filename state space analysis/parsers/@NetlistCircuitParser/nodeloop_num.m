% function [A,B,C,D,I,HtempAB,dependsAB,HtempCD,savedCD,StateNamesAB,StateNamesCD,OutputNames,DependentNames,SortedTree1,SortedCoTree1,ConstantNames,OrderedNamesnum,almost_H] = nodeloop_num(obj,NL,NLnets)
function [almost_H] = nodeloop_num(obj,NL,NLnets)
% nodeloop creates the state matrices A,B,C,D from a specific node input
% matrix
%
% Notes: Needs to be more general, broken up into a few more functions,
%
%

%% Sort numeric index values 
% This section matches the user given values to all variables from the
% topology level to allow numerical solving through a lookup table
% 
[~,IA,IB] = intersect(NLnets(:,1),obj.Component_Values(:,1));
obj.index = zeros(size(NLnets,1),1);
obj.index(IA) = IB;




 D_L = obj.Cutset;
 SortedTree = obj.SortedTree_cutloop;
 SortedCoTree = obj.SortedCoTree_cutloop;

% [row,col]=size(D_L);

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


% hey = find(SortedTree(:,1)==numR,1); % find first index of resistor (R) in tree
% if isempty(hey) % if there are no R in tree then set to zero
%     R = row+1;
%     Z_R = 0;
%     Y_R = 0;
% else % if there are resistors in tree
%     % Create symbolic values for all resistor names:
%     R = hey;
%     %syms(NLnets(SortedTree(hey:end,4),1));
%     %syms 'Z_R';
%     % Create numeric resistor matrix for Tree:
%     for i = 1:1:row-R+1
%         Z_R(i,i) = obj.Component_Values{obj.index(SortedTree(R+i-1,4)),2};
%         % Z_R(i,i) = NLnets(SortedTree(R+i-1,4),1); % What the
%         % symbolic vesion was
%     end
%     % Find inverse of symbolic resistor matrix for Tree:
%     Y_R = inv(Z_R);
% end
Rlocs = SortedTree(:,1)==numR;
Relements = obj.index(SortedTree(Rlocs,4));
Z_R = diag([obj.Component_Values{Relements,2}]);
Y_R = inv(Z_R);
if(isempty(Z_R))
    Z_R = 0;
    Y_R = 0;
end




% ya = find(SortedCoTree(:,1)==numG,1,'last'); % find last index of resistor (R) in tree
% if isempty(ya) % if there are no resistors in cotree then set to zero
%     G = 0;
%     Y_G = 0;
%     Z_G = 0;
% else % if there are resistors in cotree
%     % Create symbolic values for all resistor names:
%     G = ya;
%     % syms(NLnets(SortedCoTree(1:ya,4),1));
%     % syms 'Z_G';
%     % Create symbolic resistor matrix for CoTree:
%     for i = 1:1:G
%         Z_G(i,i) = obj.Component_Values{obj.index(SortedCoTree(i,4)),2};
%         % Z_G(i,i) = NLnets(SortedCoTree(i,4),1);
%     end
%     % Find inverse of symbolic resistor matrix for CoTree:
%     Y_G = inv(Z_G);
% end

Glocs = SortedCoTree(:,1)==numG;
Gelements = obj.index(SortedCoTree(Glocs,4));
Z_G = diag([obj.Component_Values{Gelements,2}]);
Y_G = inv(Z_G);
if(isempty(Z_G))
    Z_G = 0;
    Y_G = 0;
end

%% Parsing D_L Matrix
% D_EG = D_L(1:R-1,1:G);
% D_EJ = D_L(1:R-1,G+1:col);
% D_RG = D_L(R:row,1:G);
% D_RJ = D_L(R:row,G+1:col);
D_EG = D_L(~Rlocs,Glocs);
D_EJ = D_L(~Rlocs,~Glocs);
D_RG = D_L(Rlocs,Glocs);
D_RJ = D_L(Rlocs,~Glocs);

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
%     almost_H = [H_EE,H_EJ;H_JE,H_JJ];
elseif Flag_RJ == 1 && Flag_RG == 1 && Flag_EG == 0 && Flag_EJ == 0
    % Correct H_JJ if there are no resistors in the tree
    H_EE = -(D_EG*(Z\(D_EG.')));%-(D_EG*inv_Z*(D_EG.'));
    H_EJ = -D_EJ;
    H_JE = (D_EJ.');
    [~,H_EJ_col]=size(H_EJ);
    [H_JE_row,~]=size(H_JE);
    H_JJ = zeros(H_JE_row,H_EJ_col);
%         almost_H = [H_EE,H_EJ;H_JE,H_JJ];
elseif  Flag_EG == 1 && Flag_RG == 1 && Flag_RJ == 0 && Flag_EJ == 0
    % Correct H_EE if there are no resistors in the cotree
    H_JJ = -((D_RJ.')*inv_Y*D_RJ);
    H_EJ = -D_EJ;
    H_JE = (D_EJ.');
    [H_EJ_row,~]=size(H_EJ);
    [~,H_JE_col]=size(H_JE);
    H_EE = zeros(H_EJ_row,H_JE_col);
%             almost_H = [H_EE,H_EJ;H_JE,H_JJ];
elseif Flag_RJ == 1 && Flag_RG == 1 && Flag_EG == 1 && Flag_EJ == 0
    H_EJ = -D_EJ;
    H_JE = (D_EJ.');
    [H_EJ_row,H_EJ_col]=size(H_EJ);
    [H_JE_row,H_JE_col]=size(H_JE);
    H_JJ = zeros(H_JE_row,H_EJ_col);
    H_EE = zeros(H_EJ_row,H_JE_col);
%                 almost_H = [H_EE,H_EJ;H_JE,H_JJ];
else
    error('Unknown combination of elements')
end



%% Find H and s Matrix

% H matrix including voltage sources
almost_H = [H_EE,H_EJ;H_JE,H_JJ];
% [H,s]=obj.hybridparse(almost_H,SortedTree,SortedCoTree);
% 
% % Functions to find outputs
% [A,B,C,D,I,HtempAB,dependsAB,StateNamesAB,OutputNames,DependentNames,ConstantNames,OrderedNamesnum]=obj.loopfixAB_num(H,s,NLnets,SortedTree,SortedCoTree);
% 
% [C,D,HtempCD,savedCD,StateNamesCD]=obj.loopfixCD_num(A,B,C,D,H,s,NLnets,SortedTree,SortedCoTree);
% 
% SortedTree1=zeros(2*size(obj.NL,1),5,1);
% SortedCoTree1=zeros(2*size(obj.NL,1),5,1);
% SortedTree1(1:size(SortedTree,1),1:5) = SortedTree;
% SortedCoTree1(1:size(SortedCoTree,1),1:5) = SortedCoTree;
end

