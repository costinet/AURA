function [] = transformers(obj,mutual,value,position)
%TRANSFORMERS creates dependent sources in place of inductors given in K
%statements
%   Transfomers takes an input of the mutual

% Currently there is no support for a coupling coefficient of a
% non-one value. Leakage and magnetizing inductance must be added as
% separate inductor to the circuit by the designer to simulate

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


% position: rows are different k statements and columns are the index in the
% NL where the inductors are for a particular k statement

NL = obj.NL;
NLnets = obj.NLnets;

D = size(position);
Values = [];


% Replaces L indexes involved in k statements with either a dependent
% voltage or current source
for i = 1:1:D(1)
    for j = 1:1:D(2)
        p = position(i,j);
        if j==1
            NL(p,1) = numBI;
            NLnets(p,:) = [strcat(NLnets(p,1),'_T'), NLnets(p,2:end)];
            Values(end+1) = obj.sip2num(NLnets{p,4});

        else
            NL(p,1) = numBV;
            NLnets(p,:) = [strcat(NLnets(p,1),'_T'), NLnets(p,2:end)];
            Values(end+1) = obj.sip2num(NLnets{p,4});
        end
    end
end

%% Create K Matrix


% i REALLY  dont know what this is maybe related to the numerical
% update stuff i was working on in the fall

%{
shortcellarray = mutual;
longcellarray  = obj.Component_Values(:,1);
index  = zeros(length(longcellarray),1);
for ind =1:length(shortcellarray)
    for ind2 = 1:length(longcellarray)
        truth = strcmp(longcellarray{ind2}, shortcellarray{ind});
        if truth
            index(ind2) = ind;
            break
        end
    end
end
%}

syms(mutual);
syms 'mutual2';





for i = 1:1:length(Values)
    mutual2(i) = mutual(i);
end

% Need to find the values of transformer values


% Note on Transformers:
%  Np/Ns = Vp/Vs
%  For two legs transformer: Vp*Ip = Vs*Is
%  (Np/Ns)^2 = Lp/Ls
%  Np/Ns = sqrt(Lp/Ls)
%  nI+nI+nI=0

% First Inductor Listed is Primary Coil

for i = 2:1:length(Values)
    %Turns(i,i-1) = sqrt(Values(1)/Values(i));
    Turns_syms(i-1,length(Values)) = sqrt(mutual2(i)/mutual2(1));
    Turns_syms(length(Values),i-1) = -sqrt(mutual2(i)/mutual2(1));
end

K = Turns_syms;



obj.NL = NL;
obj.NLnets = NLnets;
obj.K = K;
J = 795723;

end % That's all Folk's
