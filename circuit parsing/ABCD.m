function [StateName] = ABCD(filename)
% ABCD creates takes a NETlist file from LTSpice and creates the
% associated ABCD matrices.

%{
% fileID = -1;
% errmsg = '';
% while fileID < 0 
%    disp(errmsg);
%    filename = input('Open file: ', 's');
%    [fileID,errmsg] = fopen(filename);in
% end
%}

%% Read in file:
[NLwhole,NLnets,NL,K]=read_file(filename);


%% Find Diodes


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


switches = [];
Diodes = [];
for i = 1:1:length(NL)
    if NL(i,1) == numD 
        switches(end+1) = i;
        Diodes(end+1,:)  = NL(i,2:3);
    end
    if NL(i,1) == numM
        switches(end+1) = i;
        Diodes(end+1,:)  = [NL(i,3),NL(i,2)];
    end
end

%{ 
'Diodes' contains the position for all of the diodes in the circuit.
These are sorted by:
 Anodes in the 1st column
 Cathodes in the 2nd column
Therefore for every row, there can not be a positive voltage from the node
in the 1st column to the node in the second column
%}

%% Get binary representation of number of states to change R and C for D and M

number_of_states = 2^length(switches);
bin=de2bi(0:number_of_states-1);
state = bin;

% state = bin+2;
% ST = [];

% for i = 1:1:number_of_states
%     for j = 1:1:length(switches)
%         NL(switches(j),1) = state(i,j);
%         % each row is a different state
%     end
%      nodeloop(NL,NLnets);
%     %[NL1(:,:,i),ST(:,:,i)] = Tree(NL,NLnets);
% 
% end

%% Cycle through all possible states

for i = 1:1:number_of_states

    [NewNL,NewNLnets]=states(NL,NLnets,state,i,switches);
    %circuitplot(NewNL,NewNLnets);
    [A(:,:,i),B(:,:,i),C(:,:,i),D(:,:,i),StateName(:,i)] = nodeloop(NewNL,NewNLnets,K);

    J = 9572839;
end
%{
    C = [eye(4)];
    D = [zeros(4,1)];
    L3 = 16e-6;
    C1 = 40e-6;
    R1 = 10;
    M1_C = 1e-9;
    D1_C = 1e-9;
    M1_R = 0.01;
    D1_R = 0.01;
    R2 = 0.001;
    
    % Switch cases for Flyback converter
    
    switch i
        case 1 % M and D off
            %  Inital condidtions
            % V_D1  V_C1 L_3  V_M1 
            
            X = [25 25 0 25];
                        
        case 2 % M on D off
            
            X = [25 25 0 25];
            
        case 3 % M off D on
            
            X = [25 25 0 25];
            
        case 4 % M and D on
            
             X = [25 50 0 -25];
            
    end
    
    
    %sys  = ss(eval(A),eval(B),eval(C),eval(D)); % Create state space
    sys  = ss(eval(A(:,:,i)),eval(B(:,:,i)),C,D);
    sys.StateName = StateName(:,i);
    sys.OutputName = StateName(:,i);
    
    % Create input and time:
    t = linspace(0,10e-6,100000);
    u = 25*ones(size(t));
    
    Y = lsim(sys,u,t,X);
    
    figure
    h = lsimplot(sys,u,t,X);
    p = getoptions(h);
    p.YLim(1) = {[min(Y(100:end,1)) max(Y(100:end,1))]};
    p.YLim(2) = {[min(Y(100:end,2)) max(Y(100:end,2))]};
    p.YLim(3) = {[min(Y(100:end,3)) max(Y(100:end,3))]};
    p.YLim(4) = {[min(Y(100:end,4)) max(Y(100:end,4))]};
    setoptions(h,p);
    
%}
    J = 5782975892;
    
end