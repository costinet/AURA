function [C,D,Htemp,saved,OutNames] = loopfixCD_num(obj,A,B,C,D,H,s,NLnets,SortedTree,SortedCoTree)
%loopfixCD Calculates C and D
%   Detailed explanation goes here

%     _   _   _  ____    _    
%    / \ | | | |/ _  |  / \   
%   / _ \| | | | (_| | / _ \  
%  / ___ | |_| |> _  |/ ___ \ 
% /_/   \_\___//_/ |_/_/   \_\


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


% Create identity matrix along with H and almost H matrix in order to find
% correct orientation of state and output matrix
[H_row,~] = size(H); % Get size of Hybrid Matrix without s
%[almostH_row,~] = size(almost_H); % Get size of full Hybrid Matrix
Htemp = [eye(H_row),H,s]; % Create Mx = Ax + Bu form
%Htemp = sym(Htemp);
%almostHtemp = [eye(almostH_row),almost_H]; % Create Mx = Ax + Bu form

temps = [SortedTree(:,:);SortedCoTree(:,:)];
cir = temps(:,1)~=numR & temps(:,1)~=numG; % find position of elements in almost H
cir_state = find(temps(:,1)~=numR & temps(:,1)~=numG & temps(:,1)~=numE & temps(:,1)~=numJ & temps(:,1)~=numEB & temps(:,1)~=numJB); % find postitin of elements in H

OrderedNamesnum = temps(cir_state,4); % Find the index for the output state names
OutputNames = NLnets(OrderedNamesnum,1); % Find the list of output names
DependentNames = {};

OrderedNameselement = temps(cir_state,1);
loop = length(OrderedNameselement)+1; % Set while loop index

% Output String Names
OutNames=NLnets(temps(((temps(:,1)==3)|(temps(:,1)==10)),4),1);


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


i = 1;
j = 0;
k = 0;
%syms 'saved';

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
        depends(j,:)=Htemp(i,H_row2+1:2*H_row2); % Index forward to get the A matrix sandwitched between the I and B matrix
        depends(j,i) = 0;

        % Htemp(1:H_row2,1:H_row2) = Htemp(1:H_row2,1:H_row2) + repmat(depends(j,:),H_row2,1).*Htemp(:,i);
        Htemp(:,i+H_row2) = [];
        if j == 1
            saved=Htemp(:,i);
        end
        saved(:,j) = Htemp(:,i);
        saved(i,:) = [];
        Htemp(i,:)=[];
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

% Delete measure states
Htemp(:,2*(H_row2-1)-sum(SortedCoTrees(:,1)==10)+1:2*(H_row2-1))=[];
Htemp(:,H_row2:H_row2+sum(SortedTrees(:,1)==3)-1)=[];

if ~isempty(C)

    if j~=0 
        % This takes the currents and voltages that were solved for
        % previously in loopfixAB.m that should be the opposite of the
        % state vector, so for C its I and for L its V. From there for
        % those state that are dependent could have a current
        % component for capacitors so add back in the current that was
        % found in the previous C and D matrix
        
        for i = 1:1:size(saved,2)
            Htemp(:,H_row2:end-size(s,2)) = Htemp(:,H_row2:end-size(s,2)) + repmat(C(end-size(saved,2)+i,1:end-size(saved,2)),H_row2-1,1).*-saved(:,i);
            Htemp(:,end+1-size(s,2):end) = Htemp(:,end+1-size(s,2):end) + repmat(D(end-size(saved,2)+i,:),H_row2-1,1).*-saved(:,i);
        end
    end
    
    % Delete state outputs:
    Htemp(sum(SortedTrees(:,1)==3)+1:H_row2-1-sum(SortedCoTrees(:,1)==10),:) = [];
    Htemp(:,sum(SortedTrees(:,1)==3)+1:H_row2-1-sum(SortedCoTrees(:,1)==10)) = [];


    Htemp = rref(Htemp);

    C = [Htemp(:,size(Htemp,1)+1:end-size(s,2)),zeros(size(Htemp,1),j)];
    D = Htemp(:,end-size(s,2)+1:end);
    
    if j == 0
        saved = [];
    end
end
end
