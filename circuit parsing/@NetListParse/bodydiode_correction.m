function [state] = bodydiode_correction(obj,switches,state)
%bodydiode_correction eliminated impossible states for body diodes and
%FETs such as the FET and the body diode being on at the same time
%   Detailed explanation goes here
%
%
% This code fixes the state matrix to correct for impossible
% states such as the FET and the body diode being on at the same time
%
% The output STATE and SWITCHES matrices is organized as follows:
%
%
% switches  = [3 3 6];
% The numbers 3 and 6 indicate the row position of switches (either
% FETs or diodes) in the netlist NL. There are two 3s because it
% indicates a FET position to represent the possibility of the FET
% being on (the first 3) and the body diode of the FET being on (the second 3)
%
%
% state =
%
%      0     0     0
%      1     0     0
%      0     1     0
%      0     0     1
%      1     0     1
%      0     1     1
%
% Not true:
% state gives all of the possible state combinations that could exist.
% This code corrects invalid states. An example is that in the above
% example there cannot be a 1 in both the 1st and 2nd column in any
% row because a FET cannot be on and have its body diode conduct at
% the same time
%
% bd_state and off_state is shows the state matrix index for all
% states being the same excpet for the columns of the indexed state
%
%
% Will later incorporate both body diode and FET being on at the same
% time due to high negative current on the FET that will exceed the
% forward voltage of the body diode
%
%     _   _   _  ____    _
%    / \ | | | |/ _  |  / \
%   / _ \| | | | (_| | / _ \
%  / ___ | |_| |> _  |/ ___ \
% /_/   \_\___//_/ |_/_/   \_\
%

% Index values for components (will pass as variable... eventually):
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


NL=obj.NL;
n = size(state,1);
k  = 0;

% This takes out the state where the body diode and FET are ON
% This will later be taken out because this state is possible****
while k<n
    k = k+1;
    for j = 1:1:size(state,2)-1
        if state(k,j)==1 && state(k,j+1)==1 && switches(j)==switches(j+1)
            state(k,:) = [];
            n = n-1;
            k = k-1;
            break
        end
    end
end

bd_state2 = [zeros(size(state))]; % Initialize matrix
bd_state = [zeros(size(state))]; % Initialize matrix
off_state = [zeros(size(state))]; % Initialize matrix

% These nested for loops find the associated diode on and off states for each state combination solved for FETs
for j = 1:1:length(switches)-1
    for i = 1:1:size(state,1)
        if switches(j)==switches(j+1) && (NL(switches(j),1)==numM) % IF the switch numbers are the same (its a FET) and NL says its a FET
            
            state_test = state(i,:);
            state_test(j:j+1) = [0 1]; % For FET selected want to find where body diode is ON and FET is OFF (this will need to later be adjusted if body diode and FET are on at the same time)
            state_test2 = state(i,:);
            state_test2(j:j+1) = [0 0]; % For FET selected want to find where body diode is OFF and FET is OFF (this will need to later be adjusted if body diode and FET are on at the same time)
            
            
            bd_state(i,j) = find(sum(repmat(state_test,size(state,1),1)==state,2)==size(state,2)==1);
            bd_state(i,j+1) = bd_state(i,j);
            
            off_state(i,j) = find(sum(repmat(state_test2,size(state,1),1)==state,2)==size(state,2)==1);
            off_state(i,j+1) = off_state(i,j);
            
            
            % Cycle through and find correct states
            %{
            for k = 1:1:size(state,1)
                if state_test == state(k,:)
                    bd_state(i,j) = k;
                    bd_state(i,j+1) = k;
                end
                if state_test2 == state(k,:)
                    off_state(i,j) = k;
                    off_state(i,j+1) = k;
                end
            end
            %}
        end
        
        
        
        % These nested for loops find the associated diode on and off states for each state combination solved for diodes
        if (NL(switches(j),1)==numD)
            
            state_test = state(i,:);
            state_test(j) = [1];
            state_test2 = state(i,:);
            state_test2(j) = [0];
            
            bd_state(i,j) = find(sum(repmat(state_test,size(state,1),1)==state,2)==size(state,2)==1);
            
            off_state(i,j) = find(sum(repmat(state_test2,size(state,1),1)==state,2)==size(state,2)==1);
            
            %{
            for k = 1:1:size(state,1)
                if state_test == state(k,:)
                    bd_state(i,j) = k;
                end
                if state_test2 == state(k,:)
                    off_state(i,j) = k;
                end
            end
            %}
        end
    end
end


% This is an exception due to the method of indexing used in the above for loop
% If the last switch in SWITCHES is a diode then was skipped in the previous step and must be included here
if (NL(switches(end),1)==numD)
    j = length(switches);
    for i = 1:1:size(state,1)
        state_test = state(i,:);
        state_test(j) = [1];
        state_test2 = state(i,:);
        state_test2(j) = [0];
        
        bd_state(i,j) = find(sum(repmat(state_test,size(state,1),1)==state,2)==size(state,2)==1);
        
        off_state(i,j) = find(sum(repmat(state_test2,size(state,1),1)==state,2)==size(state,2)==1);
        
        %{
        for k = 1:1:size(state,1)
            if state_test == state(k,:)
                bd_state(i,j) = k;
            end
            if state_test2 == state(k,:)
                off_state(i,j) = k;
            end
        end
        %}
    end
    
end

obj.BD_state = bd_state;
obj.BD_OFF_state = off_state;

end % That's all Folks
