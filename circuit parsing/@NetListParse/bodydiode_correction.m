function [state] = bodydiode_correction(obj,switches,state)
%bodydiode_correction eliminated impossible states for body diodes and
%FETs such as the FET and the body diode being on at the same time
%   Detailed explanation goes here
%
%
% This code fixes the the state matrix to correct for impossible
% states such such as the FET and the body diode being on at the same time
%
% The ouput STATE and SWITCHES matricies is orgainized as follows:
%
%
% swithces  = [3 3 6];
% The numbers 3 and 6 indicate the row position of swithces (either
% FETs or diodes) in the netlist NL. There are two 3s becuase it
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
% state gives all of the possible state combinations that could exist.
% This code corrects invalid states. An example is that in the above
% example there cannot be a 1 in both the 1st and 2nd colomun in any
% row becuase a FET cannot be on and have its body diode conduct at
% the same time
%
%
%
%     %%%%%%   %      %  %%%%%%%    %%%%%%
%    %      %  %      %  %      %  %      %
%    %      %  %      %  %      %  %      %
%    %%%%%%%%  %      %  %%%%%%%   %%%%%%%%
%    %      %  %      %  %%        %      %
%    %      %  %      %  % %       %      %
%    %      %  %      %  %  %      %      %
%    %      %  %      %  %   %     %      %
%    %      %   %    %   %    %    %      %
%    %      %    %%%%    %     %   %      %
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


bd_state = [zeros(size(state))]; % Initallize matrix



for j = 1:1:length(switches)-1
    if switches(j)==switches(j+1) && (NL(switches(j),1)==numM)
        for i = 1:1:size(state,1)
            state_test = state(i,:);
            state_test(j:j+1) = [0 1];
            for k = 1:1:size(state,1)
                if state_test == state(k,:)
                    bd_state(i,j) = k;
                    bd_state(i,j+1) = k;
                    break
                end
            end
        end
        j = j+1; % On purpose to speed up (No need to check for body diode since ON state is the flag for FETs)
    end
    if (NL(switches(j),1)==numD)
        for i = 1:1:size(state,1)
            state_test = state(i,:);
            state_test(j) = [1];
            for k = 1:1:size(state,1)
                if state_test == state(k,:)
                    bd_state(i,j) = k;
                end
            end
        end
    end
end

end % That's all Folks

