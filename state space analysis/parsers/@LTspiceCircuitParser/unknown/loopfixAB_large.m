function [A,B,C,D] = loopfixAB_large(obj,Htemp,depends,OutputNames,DependentNames)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Find the number of dependent states
j = size(DependentNames,1);

% Find the number of states (used for indexing)
[H_row2,~] = size(Htemp);
H_row2 = H_row2+1;

% RREF Htemp with custom tolerance value
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

for i = 1:1:j
    outdependstate = depends(i,:)'.*Htemp(:,H_row2:2*(H_row2-1));
    outdependsconst = depends(i,:)'.*Htemp(:,(2*(H_row2-1))+1:end);
    outstatedependsconst(i,:) = sum(outdependsconst).*DependentNames(i);
    outstatedepends(i,:)=sum(outdependstate).*DependentNames(i);
    % add all columns of matrix to get equation that goes in state equation
end

% Create C and D matrix
C = [OutputHtemp(:,H_row2:2*(H_row2-1)).*OutputNames,zeros(H_row2-1,j);outstatedepends(:,:),zeros(j,j)];
D = [OutputHtemp(:,(2*(H_row2-1))+1:end).*OutputNames;outstatedependsconst];


end % That's all Folks
