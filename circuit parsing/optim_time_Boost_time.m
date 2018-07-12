function [] = optim_time_Boost_time(parse)
%optim_time_Boost tests the time it takes for eval() to update the
%number of the listed elements values given in changed_values
% It takes a NetListParse class that has been solved using Test_Parse.m

% If you want to have a record of all lines of code run then set to true
% and check log file testeval.txt
log_output = false;

% Enable diary and echo for log
if log_output
    delete testeval_time.txt
    diary testeval_time.txt
    echo on all
end

% Set the number of iterations to update
number_of_iterations = 1;

% Variable ranges for Boost converter example
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

% Set random variables based on ranges
L1x = (L1b-L1a)*rand(number_of_iterations,1)+L1a;
C1x = (C1b-C1a)*rand(number_of_iterations,1)+C1a;
R1x = (R1b-R1a)*rand(number_of_iterations,1)+R1a;
M1_Cx = (M1_Cb-M1_Ca)*rand(number_of_iterations,1)+M1_Ca;
M1_Rx = (M1_Rb-M1_Ra)*rand(number_of_iterations,1)+M1_Ra;
D1_Cx = (D1_Cb-D1_Ca)*rand(number_of_iterations,1)+D1_Ca;
D1_Rx = (D1_Rb-D1_Ra)*rand(number_of_iterations,1)+D1_Ra;

% Initalize evaltime
evaltime = [];

% Set SortedTree and SortedCoTree from class
SortedTree = parse.SortedTree;
SortedCoTree = parse.SortedCoTree;

% Check to see if NetListParse class has been solved symbolically
if isempty(parse.Asym)
    
    % Go through the number of iterations given
    for i = 1:1:number_of_iterations
        
        % Start Timer
        tic
        
        % Update element values
        L1 = L1x(i);
        C1 = C1x(i);
        R1 = R1x(i);
        M1_C = M1_Cx(i);
        M1_R = M1_Rx(i);
        D1_C = D1_Cx(i);
        D1_R = D1_Rx(i);
        
        toc
        
        % Iterate trough each possible switching combination given
        for k = 1:1:size(parse.HtempAB,3)
            tic
            HtempAB(:,:,k) = eval(parse.HtempAB(:,:,k));
            toc
            tic
            HtempCD(:,:,k) = eval(parse.HtempCD(:,:,k));
            toc
            tic
            dependsAB(:,:,k) = eval(parse.dependsAB(:,:,k));
            toc
            tic
            savedCD(:,:,k) = eval(parse.savedCD(:,:,k));
            toc
            for j = 1:1:size(parse.DependentNames(:,k),1)
                tic
                DependentNames(j,k) = eval(parse.DependentNames{j,k});
                toc
            end
            for j = 1:1:size(parse.OutputNames(:,k),1)
                tic
                OutputNames(j,k) = eval(parse.OutputNames{j,k});
                toc
            end
        end
        tic
        % Solve ABCD and update NetListParse class
        for k = 1:1:size(parse.HtempAB,3)
            [A,B,C,D] = parse.loopfixAB_large(HtempAB(:,:,k),dependsAB(:,:,k),OutputNames(:,k),DependentNames(:,k));
            [C,D] = parse.loopfixCD_large(B,C,D,HtempCD(:,:,k),savedCD(:,:,k),DependentNames(:,k),SortedTree(:,:,k),SortedCoTree(:,:,k));
            parse.Anum(:,:,k)=A;
            parse.Bnum(:,:,k)=B;
            parse.Cnum(:,:,k)=C;
            parse.Dnum(:,:,k)=D;
        end
        
        % Record time it took to go through one update
        toc
        
    end
else % If the symbolic matrix has been solve and is in class
    
    % Iterate trough each possible switching combination given
    for i = 1:1:number_of_iterations
        
        % Start Timer
        tic
        
        % Update element values
        L1 = L1x(i);
        C1 = C1x(i);
        R1 = R1x(i);
        M1_C = M1_Cx(i);
        M1_R = M1_Rx(i);
        D1_C = D1_Cx(i);
        D1_R = D1_Rx(i);
        
        toc
        
        % Update NetListParse
        for k = 1:1:size(parse.Asym,3)
            tic
            parse.Anum(:,:,k) = eval(parse.Asym(:,:,k));
            toc
            tic
            parse.Bnum(:,:,k) = eval(parse.Bsym(:,:,k));
            toc
            tic
            parse.Cnum(:,:,k) = eval(parse.Csym(:,:,k));
            toc
            tic
            parse.Dnum(:,:,k) = eval(parse.Dsym(:,:,k));
            toc

            
        end
    end
    % Record time it took to go through one update
    %evaltime(end+1) = toc;
end

% Calculate metrics for time
%

% Turn off echo and diary for log
if log_output
    echo off all
    diary off
end

end % That's all folks
