function [NL,NLnets] = states(NL,NLnets,state,i,switches) 
% this function finds the states and sets the net list for nodeloop
%

% Index values for components (will pass as variable):
numV = 1;
% BV = 2;
numC = 2;
numR = 3;
numL = 4;
numI = 5;
% BI = 6;
numD = 6;
numM = 7;
numB = 8;


for j = 1:1:length(switches)
    
    % OPEN switch is replaced with cap
    if state(i,j) == 0
        NL(switches(j),1) = numC;
        NLnets(switches(j),:) = [strcat(NLnets(switches(j),1),'_C'), NLnets(switches(j),2:end)];
    end
    
    % CLOSED switch is replaced with cap and resistor
    if state(i,j) == 1
        NL(switches(j),1) = numC;
        [row_NL,~]=size(NL);
        NL(row_NL+1,:) = [numR,NL(switches(j),2:3),row_NL+1];
        NLnets(row_NL+1,:) = [strcat(NLnets(switches(j),1),'_R'), NLnets(switches(j),2:end)];
        NLnets(switches(j),:) = [strcat(NLnets(switches(j),1),'_C'), NLnets(switches(j),2:end)];
    end
    
    if state(i,j) ~= 1 && state(i,j) ~= 0
        error('Invalid State Detected')
    end
    
end