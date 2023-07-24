function [A,B,C,D,Htemp,depends,StateNames,OutputNames,DependentNames,ConstantNames,OrderedNamesnum] = loopfixAB(obj,H,s,NLnets,SortedTree,SortedCoTree)
%loopfixAB computes the ABCD matrix
%   This function either computes symbolic Htemp matrix or the ABCD matrix
%   The

% if the number of state variables is greater than sym_comput, sym AB
% matrix will not be computed
sym_comput = 99;  %% 99 or 1

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

% Eliminate Elements that no longer exist in the Hybrid Matrix from
% the Tree and CoTree Matrices
SortedCoTrees(SortedCoTrees(:,1)==numG,:)=[];
SortedTrees(SortedTrees(:,1)==numR,:)=[];
SortedCoTrees(SortedCoTrees(:,1)==numJ,:)=[];
SortedTrees(SortedTrees(:,1)==numE,:)=[];
SortedCoTrees(SortedCoTrees(:,1)==numJB,:)=[];
SortedTrees(SortedTrees(:,1)==numEB,:)=[];

% Only include the state values for the H and s matrix:
H = H(sum(SortedTrees(:,1)==3)+1:size(SortedTrees,1)+size(SortedCoTrees,1)-sum(SortedCoTrees(:,1)==10),sum(SortedTrees(:,1)==3)+1:size(SortedTrees,1)+size(SortedCoTrees,1)-sum(SortedCoTrees(:,1)==10));
s = s(sum(SortedTrees(:,1)==3)+1:size(SortedTrees,1)+size(SortedCoTrees,1)-sum(SortedCoTrees(:,1)==10),:);

% Creat identity matrix along with H and almost H matrix in order to find
% correct orentation of state and output matrix
[H_row,~] = size(H); % Get size of Hybrid Matrix without s
Htemp = [eye(H_row),H,s]; % Create Mx = Ax + Bu form
Htemp = sym(Htemp); % Ensure Htemp is a symbolic matrix

temps = [SortedTree(:,:);SortedCoTree(:,:)];
cir = temps(:,1)~=numR & temps(:,1)~=numG; % find position of elements in almost H
cir_state = find(temps(:,1)~=numR & temps(:,1)~=numG & temps(:,1)~=numE & temps(:,1)~=numJ & temps(:,1)~=numEB & temps(:,1)~=numJB & temps(:,1)~=numEM & temps(:,1)~=numJM); % find postitin of elements in H

OrderedNamesnum = temps(cir_state,4); % Find the index for the output state names
OutputNames = NLnets(OrderedNamesnum,1); % Find the list of output names
DependentNames = {}; % Initallize dependent Names
ConstantNames = NLnets(temps((temps(:,1)==numE | temps(:,1)==numJ),4),1); % Find Input names

OrderedNameselement = temps(cir_state,1);
loop = length(OrderedNameselement)+1; % Set while loop index

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
% If the number exits in the M matrix then the entire row gets
% divided by that elements either L or C
% The position for that element however does not since it forms either
% di/dt or dv/dt

for i = 1:1:length(OrderedNameselement)
    for j = 1:1:length(OrderedNameselement)
        if Htemp(i,j)~=0
            Htemp(i,:)=Htemp(i,:)./OutputNames(j);
            Htemp(i,j) = Htemp(i,j).*OutputNames(j);
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
        
        % Correct the numerical reference from the state variables to net list 
        OrderedNamesnum(end+1) = OrderedNamesnum(i); 
        OrderedNamesnum(i) = [];
        
        j = j+1; % J is number of dependent elements
        
        % Set up so each depends row is a linear combination of the
        % independent state remaining:
        % depends does not include DC sources because dv/dt of DC source is
        % zero
        depends(j,:)=Htemp(i,H_row2+1:2*H_row2);
        depends(j,i) = 0;
        
        % Any d /dt dependence based on this states gets added back to
        % the M Matrix
        Htemp(1:H_row2,1:H_row2) = Htemp(1:H_row2,1:H_row2) + repmat(depends(j,:),H_row2,1).*repmat(Htemp(:,i),1,H_row2);
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
end



% Need to format output how i want it
[H_row2,~] = size(Htemp);
H_row2 = H_row2+1;



%% For small number of state variables solve ABCD

A  = [];
B  = [];
C  = [];
D  = [];
StateNames = [OutputNames;DependentNames];

if j==0
    depends = [];
    DependentNames = cell.empty(1,0);
end

if H_row2 < sym_comput
    
   % COMPEL_2020_AURA_TEST_DAB
    
    
    Htemp = rref(Htemp);
    
    OutputHtemp = Htemp;
    
    if j~=0
        for i = 1:1:j
            dependstate = depends(i,:)'.*Htemp(:,H_row2:2*(H_row2-1));
            dependsconst = depends(i,:)'.*Htemp(:,(2*(H_row2-1))+1:end);
            statedependsconst(i,:) = sum (dependsconst);
            statedepends(i,:)=sum(dependstate);
            % add all columns of matrix to get equation that goes in state equation
        end
        
        A = [Htemp(:,H_row2:2*(H_row2-1)),zeros(H_row2-1,j);statedepends(:,:),zeros(j,j)];
        
        B = [Htemp(:,(2*(H_row2-1))+1:end);statedependsconst];
        
        
        
        %OutputHtemp = rref(OutputHtemp);
        
        
        for i = 1:1:j
            outdependstate = depends(i,:)'.*Htemp(:,H_row2:2*(H_row2-1));
            outdependsconst = depends(i,:)'.*Htemp(:,(2*(H_row2-1))+1:end);
            outstatedependsconst(i,:) = sum(outdependsconst).*DependentNames(i);
            outstatedepends(i,:)=sum(outdependstate).*DependentNames(i);
            % add all columns of matrix to get equation that goes in state equation
        end
        
        C = [Htemp(:,H_row2:2*(H_row2-1)).*OutputNames,zeros(H_row2-1,j);outstatedepends(:,:),zeros(j,j)];
        
        D = [Htemp(:,(2*(H_row2-1))+1:end).*OutputNames;outstatedependsconst];
        
    else
        depends = [];
        A = [Htemp(:,H_row2:2*(H_row2-1))];
        B = [Htemp(:,(2*(H_row2-1))+1:end)];
        C = [Htemp(:,H_row2:2*(H_row2-1)).*OutputNames];
        D = [Htemp(:,(2*(H_row2-1))+1:end).*OutputNames];
        
    end
    
end

end % That's all Folks

