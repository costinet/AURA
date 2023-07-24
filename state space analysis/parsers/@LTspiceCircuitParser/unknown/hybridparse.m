function [H,s] = hybridparse(obj,preH,SortedTree,SortedCoTree)
%hybridparse parses the hybrid matrix
%   hybridparse takes an input of the Hybrid Matrix, Sorted Tree and Sorted
%   CoTree variables and outputs the H and s matrix for the converter
%
% Currently the states are in the following order (excluding R and G)
% __________Tree__________  _________CoTree__________
% E  E-B  E-M  E-C  E-L  R  G   J-C  J-L J-M  J-B  J
% 1   2    3    4    5   6  7    8    9   10   11  12
%
%
% Matrix (preH) at the end of this function takes the form:
%  _   _       _                      _    _   _
% | W_1 |     | H_11  H_12  H_13  H_14 |  | X_1 |  State Variables
% |     |     |                        |  |     |
% | W_2 |     | H_21  H_22  H_23  H_24 |  | X_2 |  Dependent States (Transformers)
% |     |  =  |                        |  |     |
% | W_3 |     | H_31  H_32  H_33  H_34 |  | X_3 |  Measurement Sources
% |     |     |                        |  |     |
% | W_4 |     | H_41  H_42  H_43  H_44 |  | X_4 |  Independent Sources
%  _   _       _                      _    _   _
%
% From which H and s can be solved for.
%
% For now assume measurement and dependent states are the same
% i.e. row 2 is row 3, W_2 = W_3, X_2 = X_3

% We assume that measure and dependent variables are the same. This can be
% done iff there are no dependent sources in the circuit other than
% those created by the transformer model

%     _   _   _  ____    _    
%    / \ | | | |/ _  |  / \   
%   / _ \| | | | (_| | / _ \  
%  / ___ | |_| |> _  |/ ___ \ 
% /_/   \_\___//_/ |_/_/   \_\

%% Set up variables

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

K = obj.K; % Set K

SortedTreeR = SortedTree; % Set Sorted Tree Reference
SortedCoTreeG = SortedCoTree; % Set Sorted CoTree Reference

SortedCoTree(SortedCoTree(:,1)==numG,:)=[]; % Manipulate Sorted CoTree to exclude resistors
SortedTree(SortedTree(:,1)==numR,:)=[]; % Manipulate Sorted Tree to exclude resistors

% Find the new size of each Sorted Matrix
DT = size(SortedTree);
DCT = size(SortedCoTree);

% Find the last voltage component of independnet sources and dependent
% sources
lastE = find(SortedTree(:,1)==numE,1,'last');
lastEB = find(SortedTree(:,1)==numEB,1,'last');
% lastESC = find(SortedTree(:,1)==numESC,1,'last');


% Find the first independent current source and dependent sources
% firstJOC = find(SortedCoTree(:,1)==numJOC,1,'first');
firstJB = find(SortedCoTree(:,1)==numJB,1,'first');
firstJ = find(SortedCoTree(:,1)==numJ,1,'first');

test=preH;
preH = test';

%% Sorts Indexes in preH

% Sorts Indexes in tree columns

% The idea with this code is to sort into preH by sorting based on
% taking the transverse of preH then sorting the columns, then taking
% the transverse of preH and sorthing the columns again so that the
% rows and columns are sorted correctly

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

% Sorts Indexes in cotree columns

Jcol = preH(:,firstJ+DT(1):end);

if isempty(Jcol)
JBcol = preH(:,firstJB+DT(1):end);
else
JBcol = preH(:,firstJB+DT(1):firstJ+DT(1)-1);
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

% Collect sub-matrices
const = [Ecol,Jcol];
depent = [EBcol,JBcol];
measure = [EBcol,JBcol];
state = [E_state,J_state];

preH = [state,depent,const]';

% Sorts Indexes in tree rows

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

% Sorts Indexes in cotree rows

Jcol = preH(:,firstJ+DT(1):end);

if isempty(Jcol)
JBcol = preH(:,firstJB+DT(1):end);
else
JBcol = preH(:,firstJB+DT(1):firstJ+DT(1)-1);
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

% Collect sub-matrices
const = [Ecol,Jcol];
depent = [EBcol,JBcol];
measure = [EBcol,JBcol];
state = [E_state,J_state];


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


F = K*H_32;

if cond(eye(size(F))-F) == Inf
    error('Singular Matrix when trying to solve for transformer integration in solver.')
end


% If F is empty, meaning there are no transformers in the circuit, then H
% and s matrix do not need to be calculated.

if isempty(F)
    H = H_11;
    s = H_14;
else
    
    H = H_11+H_12*((eye(size(F))-F)^-1)*K*H_31;
    s = H_14+H_12*((eye(size(F))-F)^-1)*K*H_34;
end

end % That's all Folks
