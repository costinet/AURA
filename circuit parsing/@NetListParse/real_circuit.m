function [check] = real_circuit(obj,X,Y)
% real_circut checks to ensure there the initial condition of states is a
% valid electrical circuit

% Input values include X and Y matrix from a lsim solve using ABCD matrix
% and NetListParse class

% The code compares the initial conditions found in X (Identity matrix for
% C) and then compares them to their respective voltage or current
% measurements given in lsim solve (Y matrix)


check = true; % Flag if check passes
index = size(obj.OutputNames,1); % Number of Output Names
[dependname] = strtok(obj.DependentNames(:,1),'_'); % Get dependent element names without '_'
[measurename,remain] = strtok(obj.OutputNamesCD(:,1)); % Get measurement names and type of measurement Voltage (V) or Current (A)
rounding = 9; % Set rounding value for check

for i = 1:1:size(dependname,1) % Iterate through number of dependent states
    elementcheck = false; % Set flag to see if check was performed
    for j = 1:1:size(measurename,1) % loop through all measurements (Y matrix
        if strcmp(dependname{i},measurename{j}) % if there is a match between the dependent name and the measurement name
            if contains(dependname{i},'L') % if it is an inductor dependent state
                if strcmp(remain{j},' A') % if it is a current measurement
                    elementcheck = true; % mark that element check as performed
                    if round(Y(1,j),9)==round(X(1,index+i),rounding) % If initial condition matches computed initial condition X(1) = Y(1)
                        %fprintf('Inital Condition %s is correct\n',obj.DependentNames{i,1})
                    else
                        warning('The intial conditions given violate either KVL or KCL')
                        fprintf('State Variable: %s \n   Given: %9f\n   Implementation: %9f\n',obj.StateNames{i+index,1},X(1,index+i),Y(1,j))
                        difference = X(1,i+index)-Y(1,j);
                        fprintf('   Difference = %5e\n',difference)
                        check = false;
                    end
                end
            else % if it is a capacitive dependent state
                if strcmp(remain{j},' V') % if it is a voltage measurement
                    elementcheck = true; % mark that element check as performed
                    if round(Y(1,j),9)==round(X(1,index+i),rounding) % If inital condition matches computed initial condition X(1) = Y(1)
                        %fprintf('Initial Condition %s is correct\n',obj.DependentNames{i,1})
                    else
                        warning('The intial conditions given violate either KVL or KCL')
                        fprintf('State Variable: %s \n   Given: %9f\n   Implementation: %9f\n',obj.StateNames{i+index,1},X(1,index+i),Y(1,j))
                        difference = X(1,i+index)-Y(1,j);
                        fprintf('   Difference = %5e\n',difference)
                        check = false;
                    end
                end
            end
        end
    end
    if ~elementcheck
        fprintf('State %s was unable to be checked\n',obj.StateNames{i+index,1})
    end
end
end % That's all Folks
