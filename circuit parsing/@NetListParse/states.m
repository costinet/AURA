function [NL,NLnets] = states(obj,state,i,switches)
% states finds the states of a converter and sets the net list for nodeloop
%

% Index values for components (will pass as variable):
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

NL = obj.NL;
NLnets = obj.NLnets;
ON_States = cell(size(switches,2),1);
OFF_States = cell(size(switches,2),1);
ON_index = 1;
OFF_index = 1;


for j = 1:1:length(switches)

    % OPEN switch is replaced with cap
    if state(i,j) == 0
        NL(switches(j),1) = numC;
        NLnets(switches(j),:) = [strcat(NLnets(switches(j),1),'_C'), NLnets(switches(j),2:end)];
        OFF_States(OFF_index,1) = NLnets(switches(j),1);
        OFF_index = OFF_index+1;
    end

    % CLOSED switch is replaced with cap and resistor
    if state(i,j) == 1
        NL(switches(j),1) = numC;
        [row_NL,~]=size(NL);
        NL(row_NL+1,:) = [numR,NL(switches(j),2:3),row_NL+1];
        NLnets(row_NL+1,:) = [strcat(NLnets(switches(j),1),'_R'), NLnets(switches(j),2:end)];
        NLnets(switches(j),:) = [strcat(NLnets(switches(j),1),'_C'), NLnets(switches(j),2:end)];
        ON_States(ON_index,1) = NLnets(switches(j),1);
        ON_index = ON_index+1;
    end
    
    if state(i,j) ~= 1 && state(i,j) ~= 0
        error('Invalid State Detected')
    end

end

obj.ON_States(:,i)=ON_States;
obj.OFF_States(:,i)=OFF_States;


end % That's all Folks
