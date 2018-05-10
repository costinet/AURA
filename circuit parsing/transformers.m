function [NLnets,NL,K] = transformers(NLnets,NL,mutual,value,position)
%transformers creates dependent sources in place of inductors given in K
%statements
%   Transfomers takes an input of 

% Set index
% Name:       V BV MV  C  R  L MI BI  I  D  M
% Identifier: 1  2  3  4  5  6  7  8  9 10 11

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

% position: rows are different k statments and columns are the index in the
% NL where the inductors are for a particular k statement

D = size(position);
Values = [];

% Replaces L indexes involved in k statements with either a dependent
% voltage or current source
for i = 1:1:D(1)
    for j = 1:1:D(2)
        p = position(i,j);
        if j==1
            NL(p,1) = numBV;
            NLnets(p,:) = [strcat(NLnets(p,1),'_BV'), NLnets(p,2:end)];
            Values(end+1) = sip2num(NLnets{p,4});
        else
            NL(p,1) = numBI;
            NLnets(p,:) = [strcat(NLnets(p,1),'_BI'), NLnets(p,2:end)];
            Values(end+1) = sip2num(NLnets{p,4});
        end
    end
end

%{
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
                Values(end+1) = sip2num(NLnets{i,4});
            else
                
                NL(i,1) = numBI;
                NLnets(i,:) = [strcat(NLnets(i,1),'_BI'), NLnets(i,2:end)];
                Values(end+1) = sip2num(NLnets{i,4});
            end
            
        end
    end

end
%}

%% Create K Matrix

syms(mutual);
syms 'mutual2';


for i = 1:1:length(Values)
    mutual2(i) = mutual(i);
end


% Note on Transformers:
%  Np/Ns = Vp/Vs
%  Np/Ns = Is/Ip
%  Vp*Ip = Vs*Is
%  (Np/Ns)^2 = Lp/Ls
%  Np/Ns = sqrt(Lp/Ls)

% First Inductor Listed is Primary Coil

for i = 2:1:length(Values)
    Turns(i,i-1) = sqrt(Values(1)/Values(i));
    Turns_syms(i,i-1) = sqrt(mutual2(1)/mutual2(i));
end


Turns(1,end+1) = Turns(2,1);
Turns_syms(1,end+1) = Turns_syms(2,1);

K = Turns_syms;


J = 795723;