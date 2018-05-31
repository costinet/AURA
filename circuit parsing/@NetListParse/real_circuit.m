function [check] = real_circuit(obj,i,u,t,X0)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Under Construction
% Currently Wrong - will always pass
check = false;
sys  = ss(obj.Anum(:,:,i),obj.Bnum(:,:,i),eye(size(obj.Anum(:,:,i))),zeros(size(obj.Bnum(:,:,i)))); % Create state space
[Y,~,X] = lsim(sys,u,t,X0);

if round(Y(1,:),15)==round(X(1,:),15)
    check = true;
else
    warning('The intial conditions given violate either KVL or KCL')
    fprintf('The inital conditions for each state\n')
    for j = 1:1:size(obj.StateNames,1)
    fprintf('State: %s \n   Given: %9f\n   Implementation: %9f\n',obj.StateNames{j,i},X(1,j),Y(1,j))
    difference = X(1,j)-Y(1,j);
    fprintf('   Difference = %5e\n',difference)
    end
end
