function [] = optim_time_Boost(parse)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% make get function function for state names of converter
% make push function for state names of converter
% function to eval new converter and solve numerically

% L1 = 16e-6;
% C1 = 40e-6;
% R1 = 10;
% M1_C = 1e-9;
% D1_C = 1e-9;
% M1_R = 0.01;
% D1_R = 0.01;

L1a = 13e-6;
L1b = 19e-6;
C1a = 35e-6;
C1b = 45e-6;
R1a = 8;
R1b = 12;
M1_Ca = 0.5e-9;
M1_Cb = 3e-9;
M1_Ra = 0.005;
M1_Rb = 0.1;
D1_Ca = 0.5e-9;
D1_Cb = 3e-9;
D1_Ra = 0.005;
D1_Rb = 0.1;

L1x = (L1b-L1a)*rand(100,1)+L1a;
C1x = (C1b-C1a)*rand(100,1)+C1a;
R1x = (R1b-R1a)*rand(100,1)+R1a;
M1_Cx = (M1_Cb-M1_Ca)*rand(100,1)+M1_Ca;
M1_Rx = (M1_Rb-M1_Ra)*rand(100,1)+M1_Ra;
D1_Cx = (D1_Cb-D1_Ca)*rand(100,1)+D1_Ca;
D1_Rx = (D1_Rb-D1_Ra)*rand(100,1)+D1_Ra;

evaltime = [];

A = parse.Asym;
B = parse.Bsym;
C = parse.Csym;
D = parse.Dsym;

SortedTree = parse.SortedTree;
SortedCoTree = parse.SortedCoTree;

for i = 1:1:100
tic
L1 = L1x(i);
C1 = C1x(i);
R1 = R1x(i);
M1_C = M1_Cx(i);
M1_R = M1_Rx(i);
D1_C = D1_Cx(i);
D1_R = D1_Rx(i);

if isempty(parse.Asym)
    for k = 1:1:size(parse.HtempAB,3)
        HtempAB(:,:,k) = eval(parse.HtempAB(:,:,k));
        HtempCD(:,:,k) = eval(parse.HtempCD(:,:,k));
        dependsAB(:,:,k) = eval(parse.dependsAB(:,:,k));
        savedCD(:,:,k) = eval(parse.savedCD(:,:,k));
        for j = 1:1:size(parse.DependentNames(:,k),1)
            DependentNames(j,k) = eval(parse.DependentNames{j,k});
        end
        for j = 1:1:size(parse.OutputNames(:,k),1)
            OutputNames(j,k) = eval(parse.OutputNames{j,k});
        end
    end
    for k = 1:1:size(parse.HtempAB,3)
        [A,B,C,D] = parse.loopfixAB_large(HtempAB(:,:,k),dependsAB(:,:,k),OutputNames(:,k),DependentNames(:,k));
        [C,D] = parse.loopfixCD_large(B,C,D,HtempCD(:,:,k),savedCD(:,:,k),DependentNames(:,k),SortedTree(:,:,k),SortedCoTree(:,:,k));
        parse.Anum(:,:,k)=A;
        parse.Bnum(:,:,k)=B;
        parse.Cnum(:,:,k)=C;
        parse.Dnum(:,:,k)=D;
    end
    
else
    for k = 1:1:size(A,3)
        parse.Anum(:,:,k) = eval(A(:,:,k));
        parse.Bnum(:,:,k) = eval(B(:,:,k));
        parse.Cnum(:,:,k) = eval(C(:,:,k));
        parse.Dnum(:,:,k) = eval(D(:,:,k));
    end
end
evaltime(end+1) = toc;
end
evaltime_mean=mean(evaltime)
evaltime_range=range(evaltime)
evaltime_max=max(evaltime)
evaltime_min=min(evaltime)
evaltime_std=std(evaltime)

end



