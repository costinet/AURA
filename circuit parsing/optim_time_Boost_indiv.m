function [] = optim_time_Boost_indiv(parse)
%optim_time_Boost_indiv tests the time it takes for eval() to update the
%number of the listed elements values given in changed_values
% It takes a NetListParse class that has been solved using Test_Parse.m

% If you want to have a record of all lines of code run then set to true
% and check log file testeval_indv.txt
log_output = false;

% Enable diary and echo for log
if log_output
    delete testeval_indv.txt
    diary testeval_indv.txt
    echo on all
end

% Set the component values that are updated at each iteration
changed_values = {'C1','L1'};
number_of_iterations = 100;

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

% Initialize element values
L1 = L1x(1);
C1 = C1x(1);
R1 = R1x(1);
M1_C = M1_Cx(1);
M1_R = M1_Rx(1);
D1_C = D1_Cx(1);
D1_R = D1_Rx(1);

% Check to see if NetListParse class has been solved symbolically
if isempty(parse.Asym)

    % Set up initial values for numerical solve
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


    % Set up flags to determine if a matrix is updated based on values
    % given in changed_values
    flagHtempAB=sum(ismember(arrayfun(@char, symvar(parse.HtempAB), 'uniform', 0),changed_values));
    flagHtempCD=sum(ismember(arrayfun(@char, symvar(parse.HtempCD), 'uniform', 0),changed_values));
    flagdependsAB=sum(ismember(arrayfun(@char, symvar(parse.dependsAB), 'uniform', 0),changed_values));
    flagsavedCD=sum(ismember(arrayfun(@char, symvar(parse.savedCD), 'uniform', 0),changed_values));
    flagDependentNames=sum(sum(ismember(parse.DependentNames,changed_values)));
    flagOutputNames=sum(sum(ismember(parse.OutputNames,changed_values)));

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

        % Iterate trough each possible switching combination given
        for k = 1:1:size(parse.HtempAB,3)

            % If a flag is set then eval each matrix needed to solve ABCD
            if flagHtempAB>0
                HtempAB(:,:,k) = eval(parse.HtempAB(:,:,k));
            end
            if flagHtempCD>0
                HtempCD(:,:,k) = eval(parse.HtempCD(:,:,k));
            end
            if flagdependsAB>0
                dependsAB(:,:,k) = eval(parse.dependsAB(:,:,k));
            end
            if flagsavedCD>0
                savedCD(:,:,k) = eval(parse.savedCD(:,:,k));
            end
            if flagDependentNames>0
                for j = 1:1:size(parse.DependentNames(:,k),1)
                    DependentNames(j,k) = eval(parse.DependentNames{j,k});
                end
            end
            if flagOutputNames>0
                for j = 1:1:size(parse.OutputNames(:,k),1)
                    OutputNames(j,k) = eval(parse.OutputNames{j,k});
                end
            end
        end

        % Solve ABCD and update NetListParse class
        for k = 1:1:size(parse.HtempAB,3)
            [A,B,C,D] = parse.loopfixAB_large(HtempAB(:,:,k),dependsAB(:,:,k),OutputNames(:,k),DependentNames(:,k));
            [C,D] = parse.loopfixCD_large(B,C,D,HtempCD(:,:,k),savedCD(:,:,k),DependentNames(:,k),SortedTree(:,:,k),SortedCoTree(:,:,k));
            parse.Anum(:,:,k)=A;
            parse.Bnum(:,:,k)=B;
            parse.Cnum(:,:,k)=C;
            parse.Dnum(:,:,k)=D;
        end
    end
    % Record time it took to go through one update
    evaltime(end+1) = toc;

else % If the symbolic matrix has been solve and is in class

    % Set up initial values for numerical solve
    for k = 1:1:size(parse.Asym,3)
        parse.Anum(:,:,k) = eval(parse.Asym(:,:,k));
        parse.Bnum(:,:,k) = eval(parse.Bsym(:,:,k));
        parse.Cnum(:,:,k) = eval(parse.Csym(:,:,k));
        parse.Dnum(:,:,k) = eval(parse.Dsym(:,:,k));
    end

    % Set up flags to determine if a matrix is updated based on values
    % given in changed_values
    flagA=sum(ismember(arrayfun(@char, symvar(parse.Asym), 'uniform', 0),changed_values));
    flagB=sum(ismember(arrayfun(@char, symvar(parse.Bsym), 'uniform', 0),changed_values));
    flagC=sum(ismember(arrayfun(@char, symvar(parse.Csym), 'uniform', 0),changed_values));
    flagD=sum(ismember(arrayfun(@char, symvar(parse.Dsym), 'uniform', 0),changed_values));

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

        % Update NetListParse
        for k = 1:1:size(parse.Asym,3)
            if flagA>0
                parse.Anum(:,:,k) = eval(parse.Asym(:,:,k));
            end
            if flagB>0
                parse.Bnum(:,:,k) = eval(parse.Bsym(:,:,k));
            end
            if flagC>0
                parse.Cnum(:,:,k) = eval(parse.Csym(:,:,k));
            end
            if flagD>0
                parse.Dnum(:,:,k) = eval(parse.Dsym(:,:,k));
            end
        end
        % Record time it took to go through one update
        evaltime(end+1) = toc;
    end
end

% Calculate metrics for time
evaltime_mean=mean(evaltime)
evaltime_range=range(evaltime)
evaltime_max=max(evaltime)
evaltime_min=min(evaltime)
evaltime_std=std(evaltime)

% Turn off echo and diary for log
if log_output
    echo off all
    diary off
end

end % That's all Folks
