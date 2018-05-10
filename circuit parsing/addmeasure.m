function [NL,NLnets] = addmeasure(NL,NLnets)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% List of Voltages wanted to measure

% Example:
% Measure Voltage from Node 1 to Node 2 and Node 4 to Node 3 would be
% written as:
% 
% Voltage = [1 2
%            4 3];

numV = 1;
numBV = 2;
numMV = 3;
numC = 4;
numR = 5;
numL = 6;
numMI = 7;
numBI = 8;
numI = 9;
numD = 10;
numM = 11;





% if max(max(Voltage))>max(max(NL(:,2:3)))
%     toohigh = max(max(Voltage));
%     high = num2str(toohigh,'%.0f\n');
%     error(['Error: ',high,' is not a valid node number'])
% end

Voltage = {'V1'
    'L1'
    'M1'
    'D1'
    'C1'
    'R1'};


Current = {'V1'
    'L1'
    'M1'
    'D1'
    'C1'
    'R1'};

%Voltage = {'V1','M1','L1'};

% List of Currentes to Measusre
% Currents listed as going from Node 1 to Node 2 and Node 4 to Node 3
% should be listed as:
% 
% Current = [1 2
%     4 3];






% if max(max(Current))>max(max(NL(:,2:3)))
%     toohigh = max(max(Current));
%     high = num2str(toohigh,'%.0f\n');
%     error(['Error: ',high,' is not a valid node number'])
% end

% Current = {};
 Voltage = {};


%{
D = size(Voltage);

Flag = zeros(size(NL,1),1);

for i = 1:1:D(1)

Match = NL(:,2:3) == Voltage(i,:);
Matches = Match(:,1).*Match(:,2);

Flag(Matches==1,:) = true;

end

% There is a long ways to go with this one
%}

    
    
for i = 1:1:length(Voltage)
    
    row = strcmp(NLnets(:,1),Voltage(i));
    NLnets(end+1,:) = [strcat(Voltage(i),' V'), NLnets(sum(row.*NL(:,4)),2:end)];
    NL(end+1,:)= [numMI,sum(row.*NL(:,2:3)),NL(end,4)+1];
    
end


for i = 1:1:length(Current)
    c_new = strcat('New_Meas_no_',num2str(i));
    row = strcmp(NLnets(:,1),Current(i));
    new_node = max(max(NL(:,2:3)))+1;
    node = NL(sum(row.*NL(:,4)),2);
    
    New_Node = strcat(NLnets(sum(row.*NL(:,4)),2),'_Meas');
    Node = NLnets(sum(row.*NL(:,4)),2);
    
    NLnets(sum(row.*NL(:,4)),2) = New_Node;
    NLnets(end+1,:) = [strcat(Current(i),' A'),Node,New_Node,NLnets(sum(row.*NL(:,4)),4:end)];
   
    NL(sum(row.*NL(:,4)),2)=new_node;
    NL(end+1,:)= [numMV,node,new_node,NL(end,4)+1];
    
    
    
    %NLnets(end+1,:) = [strcat(Current(i),'_Current'), NLnets(sum(row.*NL(:,4)),2:end)];

    
end

%{
D=size(Voltage);
i = 1;
while i<D(1)+1
    NL(end+1,:)= [numMV,Voltage(i,:),NL(end,4)+1];
    i = i+1;
    
end

D=size(Current);
i = 1;
while i<D(1)+1
    NL(end+1,:)= [numMI,Current(i,:),NL(end,4)+1];
    i = i+1;
    
end

%}


J = 85932045;
end

