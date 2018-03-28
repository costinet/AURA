function [H,s] = hybridparse(preH,k,SortedTree,SortedCoTree)
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
numEC = 3;
numEL = 4;
numR = 5;
numG = 6;
numJC = 7;
numJL = 8;
numJB = 9;
numJ = 10;


lastE = find(SortedTree(:,1)==numE,1,'last');
lastEB = find(SortedTree(:,1)==numEB,1,'last');
% lastESC = find(SortedTree(:,1)==numESC,1,'last');



% firstJOC = find(SortedCoTree(:,1)==numJOC,1,'first');
firstJB = find(SortedCoTree(:,1)==numJB,1,'first');
firstJ = find(SortedCoTree(:,1)==numJ,1,'first');


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



H = H11+H12*((1-K*H32)^-1)*K*H31;
s = H14+H12*((1-K*H32)^-1)*K*H34;

outputArg2 = inputArg2;
end

