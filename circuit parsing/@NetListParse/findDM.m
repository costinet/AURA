function [switches] = findDM(obj)
% blah
    % blah
%% Set Variables

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

%% Find Diodes and Switches

NL = obj.NL;
NLnets = obj.NLnets;

switches = [];
Diodes = [];
for i = 1:1:size(NL,1)
    if NL(i,1) == numD
        switches(end+1) = i;
        Diodes(end+1,:)  = NL(i,2:3);
    end
    if NL(i,1) == numM
        switches(end+1) = i;
        Diodes(end+1,:)  = [NL(i,3),NL(i,2)];
    end
end

obj.Diodes = Diodes;
obj.Switches = switches;


end