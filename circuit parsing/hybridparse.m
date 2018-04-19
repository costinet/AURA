function [H,s] = hybridparse(preH,K,SortedTree,SortedCoTree)
%hybridparse parses the hybird matrix 
%   Detailed explanation goes here

%{
What it will be: Assume input tree is in the form:
% E  E-B  E-SC  E-C  E-L  R    G   J-C   J-L  J-OC   J-B   J
% 1   2    3     4    5   6    7    8     9    10    11  12

% 
% ##### What it is today: New Index Values $$$$$$$
% Branch Identification numbers:
% E  E-B  E-C  E-L  R   G   J-C  J-L  J-B  J
% 1   2    3   4    5   6    7   8     9   10


%}


% For now assume measurement and dependent states are the same rows 2 is
% row 3


% Matrix (preH) is now of the form:
%  _   _       _                      _    _   _
% | W_1 |     | H_11  H_12  H_13  H_14 |  | X_1 |  State Variables
% |     |     |                        |  |     |
% | W_2 |     | H_21  H_22  H_23  H_24 |  | X_2 |  Dependent States (Transformers)
% |     |  =  |                        |  |     |
% | W_3 |     | H_31  H_32  H_33  H_34 |  | X_3 |  Measurement Sources
% |     |     |                        |  |     |
% | W_4 |     | H_41  H_42  H_43  H_44 |  | X_4 |  Independent Sources
%  _   _       _                      _    _   _



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

SortedTreeR = SortedTree;
SortedCoTreeG = SortedCoTree;

SortedCoTree(SortedCoTree(:,1)==numG,:)=[];
SortedTree(SortedTree(:,1)==numR,:)=[];

DT = size(SortedTree);
DCT = size(SortedCoTree);

lastE = find(SortedTree(:,1)==numE,1,'last');
lastEB = find(SortedTree(:,1)==numEB,1,'last');
% lastESC = find(SortedTree(:,1)==numESC,1,'last');


% firstJOC = find(SortedCoTree(:,1)==numJOC,1,'first');
firstJB = find(SortedCoTree(:,1)==numJB,1,'first');
firstJ = find(SortedCoTree(:,1)==numJ,1,'first');

test=preH;

%% Sorts Indexes in tree

Ecol = preH(:,1:lastE);

if isempty(Ecol)
EBcol = preH(:,1:lastEB);
else
EBcol = preH(:,lastE+1:lastEB);
end

if isempty(EBcol)
    if isempty(Ecol)
        E_state = preH(:,1:DT(1));
    else
        E_state = preH(:,lastE+1:DT(1));

    end
else
    E_state = preH(:,lastEB+1:DT(1));
end

%% Sorts Indexes in cotree

Jcol = preH(:,firstJ+DT(1):end);

if isempty(Jcol)
JBcol = preH(:,firstJB+DT(1):end);
else
JBcol = preH(:,firstJB+DT(1):firstJB+DT(1)-1);
end

if isempty(JBcol)
    if isempty(Jcol)
        J_state = preH(:,DT(1)+1:end);
    else
        J_state = preH(:,DT(1)+1:DT(1)+firstJ-1);

    end
else
    J_state = preH(:,DT(1)+1:DT(1)+firstJB-1);
end

const = [Ecol,Jcol];
depent = [EBcol,JBcol];
measure = [EBcol,JBcol];
state = [E_state,J_state];

preH = [state,depent,const]';

%% Sorts Indexes in tree

Ecol = preH(:,1:lastE);

if isempty(Ecol)
EBcol = preH(:,1:lastEB);
else
EBcol = preH(:,lastE+1:lastEB);
end

if isempty(EBcol)
    if isempty(Ecol)
        E_state = preH(:,1:DT(1));
    else
        E_state = preH(:,lastE+1:DT(1));

    end
else
    E_state = preH(:,lastEB+1:DT(1));
end

%% Sorts Indexes in cotree

Jcol = preH(:,firstJ+DT(1):end);

if isempty(Jcol)
JBcol = preH(:,firstJB+DT(1):end);
else
JBcol = preH(:,firstJB+DT(1):firstJB+DT(1)-1);
end

if isempty(JBcol)
    if isempty(Jcol)
        J_state = preH(:,DT(1)+1:end);
    else
        J_state = preH(:,DT(1)+1:DT(1)+firstJ-1);

    end
else
    J_state = preH(:,DT(1)+1:DT(1)+firstJB-1);
end

const = [Ecol,Jcol];
depent = [EBcol,JBcol];
measure = [EBcol,JBcol];
state = [E_state,J_state];

preH = [state,depent,const]';



%% Organize

size_state = size(state);
size_depent = size(depent);
size_measure = size(measure);
size_const = size(const);
H_11 = state(1:size_state(2),:);
H_14 = const(1:size_state(2),:);
H_12 = depent(1:size_state(2),:);
H_31 = state(size_state(2)+1:size_state(1)-size_const(2),:); % will need to alter after measure fix
H_34 = const(size_state(2)+1:size_const(1)-size_const(2),:); % will need to alter after measure fix
H_32 = depent(size_state(2)+1:size_const(1)-size_const(2),:); % will need to alter after measure fix

%{
if isempty(lastEB) && isempty(firstJB)
    newH = [preH(:,1:lastE),preH(:,firstJ:end)];
    D=length(newH);
    H = almost_H(D(2)+1:end,D(2)+1:end);
    s = almost_H(D(2)+1:end,1:D(2));
end
    
if isempty(lastE) && isempty(firstJ)
    H_41 = 0;
    H_42 = 0;
    H_43 = 0;
    H_44 = 0;
end



if isempty(lastEB) && isempty(firstJB)
    H_21 = 0;
    H_22 = 0;
    H_23 = 0;
    H_24 = 0;
end
%}



F = K*H_32;
if isempty(F)
    H = H_11;
    s = H_14;
else
    H = H_11+H_12*((eye(size(F))-F)^-1)*K*H_31;
    s = H_14+H_12*((eye(size(F))-F)^-1)*K*H_34;
end
J = 89208923;

end

