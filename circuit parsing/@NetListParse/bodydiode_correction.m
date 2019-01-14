function [state] = bodydiode_correction(~,switches,state)
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
% FETs or diodes) in the the netlist NL. There are two 3s becuase it
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
end

