function [check] = real_circuit(obj,X,Y)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Under Construction
% Currently Wrong - will always pass
check = true;
elementcheck = false;
% sys  = ss(obj.Anum(:,:,i),obj.Bnum(:,:,i),eye(size(obj.Anum(:,:,i))),zeros(size(obj.Bnum(:,:,i)))); % Create state space
% [Y,~,X] = lsim(sys,u,t,X0);
%
% if round(Y(1,:),15)==round(X(1,:),15)
%     check = true;
%
% else
%
%     warning('The intial conditions given violate either KVL or KCL')
%     fprintf('The inital conditions for each state\n')
%     for j = 1:1:size(obj.StateNames,1)
%     fprintf('State: %s \n   Given: %9f\n   Implementation: %9f\n',obj.StateNames{j,i},X(1,j),Y(1,j))
%     difference = X(1,j)-Y(1,j);
%     fprintf('   Difference = %5e\n',difference)
%     end
%
%
%
% end

index = size(obj.OutputNames,1);


[dependname] = strtok(obj.DependentNames(:,1),'_');
[measurename,remain] = strtok(obj.OutputNamesCD(:,1));

for i = 1:1:size(dependname,1)
    
    for j = 1:1:size(measurename,1)
        
        if strcmp(dependname{i},measurename{j})
            if contains(dependname{i},'L')
                if strcmp(remain{j},' A')
                    elementcheck = true;
                    if round(Y(1,j),9)==round(X(1,index+i),9)
                        %fprintf('Inital Condition %s is correct\n',obj.DependentNames{i,1})
                    else
                        warning('The intial conditions given violate either KVL or KCL')
                        fprintf('State Variable: %s \n   Given: %9f\n   Implementation: %9f\n',obj.StateNames{i+index,1},X(1,index+i),Y(1,j))
                        difference = X(1,i+index)-Y(1,j);
                        fprintf('   Difference = %5e\n',difference)
                        check = false;
                    end
                end   
            else
                if strcmp(remain{j},' V')
                    elementcheck = true;
                    if round(Y(1,j),9)==round(X(1,index+i),9)
                        %fprintf('Inital Condition %s is correct\n',obj.DependentNames{i,1})
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

end

