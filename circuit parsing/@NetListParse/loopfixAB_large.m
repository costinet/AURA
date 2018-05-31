function [A,B,C,D] = loopfixAB_large(obj,Htemp,depends,OutputNames,DependentNames)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Need to pass H_row2 as well
% Delete ABCD
%{
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
%}

j = size(DependentNames,1);

[H_row2,~] = size(Htemp);
H_row2 = H_row2+1;

Htemp = rref(Htemp,0.9);

OutputHtemp = Htemp;


for i = 1:1:j
    dependstate = depends(i,:)'.*Htemp(:,H_row2:2*(H_row2-1));
    dependsconst = depends(i,:)'.*Htemp(:,(2*(H_row2-1))+1:end);
    statedependsconst(i,:) = sum (dependsconst);
    statedepends(i,:)=sum(dependstate);
    % add all columns of matrix to get equation that goes in state equation
end

A = [Htemp(:,H_row2:2*(H_row2-1)),zeros(H_row2-1,j);statedepends(:,:),zeros(j,j)];


B = [Htemp(:,(2*(H_row2-1))+1:end);statedependsconst];


%StateNames = [OutputNames;DependentNames];


for i = 1:1:j
    outdependstate = depends(i,:)'.*Htemp(:,H_row2:2*(H_row2-1));
    outdependsconst = depends(i,:)'.*Htemp(:,(2*(H_row2-1))+1:end);
    outstatedependsconst(i,:) = sum(outdependsconst).*DependentNames(i);
    outstatedepends(i,:)=sum(outdependstate).*DependentNames(i);
    % add all columns of matrix to get equation that goes in state equation
end

C = [OutputHtemp(:,H_row2:2*(H_row2-1)).*OutputNames,zeros(H_row2-1,j);outstatedepends(:,:),zeros(j,j)];

D = [OutputHtemp(:,(2*(H_row2-1))+1:end).*OutputNames;outstatedependsconst];


end

