function [NLnets,NL,K] = transformers(NLnets,NL,mutual,value)
%transformers creates dependent sources in place of inductors given in K
%statements
%   Detailed explanation goes here


% Need to define BV and BI numbers

numV = 1;
numBV = 2;
numC = 3;
numR = 4;
numL = 5;
numBI = 6;
numI = 7;
numD = 8;
numM = 9;

D=size(NL);

DD = size(mutual);

Values = [];

for j = 1:1:DD(2)
    for i = 1:1:D(1)
        
        emptyish = strfind(NLnets{i,1},mutual{j});
        if ~isempty(emptyish) % if there was componenets
            if j==1
                
                NL(i,1) = numBV;
                NLnets(i,:) = [strcat(NLnets(i,1),'_BV'), NLnets(i,2:end)];
            
            else
                
                NL(i,1) = numBI;
                NLnets(i,:) = [strcat(NLnets(i,1),'_BI'), NLnets(i,2:end)];

            end
           
        end
    end

end

K = 1;

J = 795723;