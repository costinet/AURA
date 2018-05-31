function [C,D] = loopfixCD_large(obj,B,C,D,Htemp,saved,DependentNames,SortedTree,SortedCoTree)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



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


SortedTrees = SortedTree;
SortedCoTrees = SortedCoTree;

SortedCoTrees(SortedCoTrees(:,1)==numG,:)=[];
SortedTrees(SortedTrees(:,1)==numR,:)=[];
SortedCoTrees(SortedCoTrees(:,1)==numJ,:)=[];
SortedTrees(SortedTrees(:,1)==numE,:)=[];
SortedCoTrees(SortedCoTrees(:,1)==numJB,:)=[];
SortedTrees(SortedTrees(:,1)==numEB,:)=[];

% Need to format output how i want it
[H_row2,~] = size(Htemp);
H_row2 = H_row2+1;

s = B; % only used for size
j = size(DependentNames,1);
for i = 1:1:size(saved,2)
    Htemp(:,H_row2:end-size(s,2)) = Htemp(:,H_row2:end-size(s,2)) + repmat(C(end-size(saved,2)+i,1:end-size(saved,2)),H_row2-1,1).*-saved(:,i);
    Htemp(:,end+1-size(s,2):end) = Htemp(:,end+1-size(s,2):end) + repmat(D(end-size(saved,2)+i,:),H_row2-1,1).*-saved(:,i);
end

% Delete state outputs:
Htemp(sum(SortedTrees(:,1)==3)+1:H_row2-1-sum(SortedCoTrees(:,1)==10),:) = [];
Htemp(:,sum(SortedTrees(:,1)==3)+1:H_row2-1-sum(SortedCoTrees(:,1)==10)) = [];



Htemp = rref(Htemp);

OutputHtemp = Htemp;


C = [Htemp(:,size(Htemp,1)+1:end-size(s,2)),zeros(size(Htemp,1),j)];
D = Htemp(:,end-size(s,2)+1:end);


end
