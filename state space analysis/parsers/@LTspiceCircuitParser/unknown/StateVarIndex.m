function [] = StateVarIndex(obj)
%STATEVARINDEX finds the index of state variables in the output vector
%   This code finds index of state variables in the output vector y = Cx+Du
%   Therefore, It find the measured voltage for capacitors and the measured
%   current for inductors.

% Disclaimer:
% Assume that all state variables have both a current and voltage
% measurement within the y matrix -- This assumption is always true as long
% as addmeasure.m is unchanged. addmeasure.m is currently set up to always
% ensure there are both current and voltage measurement elements for each
% state space.

StateNames = obj.StateNames;
OutputNames = obj.OutputNamesCD;

[measurename,remain] = strtok(OutputNames(:,1)); % Separate out the name of elements from their measured value (V or A)
[StateNames1,~] = strtok(StateNames(:,1),'_'); % Separate names from the '_C'

for i=1:1:size(StateNames1,1)
    result = strcmp(measurename,StateNames1(i));

    if contains(StateNames(i),'L')
        StateNumbers(i)=find(result==1,1,'first');
        StateNumbers_Opp(i)=find(result==1,1,'last');
    elseif contains(StateNames(i),'C')
        StateNumbers_Opp(i)=find(result==1,1,'first');
        StateNumbers(i)=find(result==1,1,'last');
    else
        fprintf('I don''t believe it.\n')
    end

end

obj.StateNumbers=StateNumbers;
obj.StateNumbers_Opposite=StateNumbers_Opp;
end % That's all Folks
