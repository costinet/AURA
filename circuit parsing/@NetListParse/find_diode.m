function [] = find_diode(obj,order)
% Function to determine where diodes are in matrix

% This code is run after ABCD in the parse class. It requires statenames
% variable created by ABCD function

% The input order is the order of states that are going to be used to test
% the converter

% Example:
% Say the converter solves the states in following order for a two
%  switch converter:
%   1) ALL OFF
%   2) FET ON
%   3) Diode ON
%   4) ALL ON
% However the time intervals for steady state are:
%   1) FET ON
%   2) ALL OFF
%   3) Diode ON
%   4) ALL OFF
% Therefore order would be:
%   order = [2 1 3 1]

% DMpos is the output
% First column if 1 if it is a Diode or FET
% Second column is 1 if Diode
% Third column is 1 if FET

% Example:
% DMpos =
%   1 0 1    % MOSFET
%   0 0 0    % Capacitor
%   0 0 0    % Inductor
%   0 1 0    % Diode

% ONOROFF is a double matrix that provides data on where switches are in
% the state matrix and what position they are in for each time interval of
% the period

% In onoroff the columns are time intervals of the period and the rows are
% the state variables

%  1 is if the FET or diode is ON
% -1 is if the FET or diode of OFF
%  0 is if the state variable is not a switch

% Example:
% Time int  |   1     2     3     4
%            ______________________
%   M1      |   1    -1    -1    -1
%   C1      |   0     0     0     0
%   L1      |   0     0     0     0
%   D1      |  -1    -1     1    -1


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



StateNames=obj.StateNames;
num=size(StateNames,1);

DMpos = zeros(num,3);

for i = 1:1:num % Iterate through number of states variables
    test=StateNames{i,1};
    if strcmp(test(1),'D') % if there is a match between the dependent name and the measurement name
        DMpos(i,:) = [1 1 0];
    elseif strcmp(test(1),'M')
        DMpos(i,:) = [1 0 1];
    else
        DMpos(i,:) = [0 0 0]; % Not necessary
        % Component is not DM
    end
end

onoroff = zeros(num,length(order));

for i = 1:1:length(order) % Iterate through number of states
    for j = 1:1:num % Iterate through number of states variables
        for k = 1:1:size(obj.ON_States,1)
            if strcmp(obj.ON_States{k,order(i)},obj.StateNames(j,1))
                onoroff(j,i)=1;
            end
            if strcmp(obj.OFF_States{k,order(i)},obj.StateNames(j,1))
                onoroff(j,i)=-1;
            end
        end
    end
end

if sum(sum(DMpos,2)==1)~=0
    fprintf('Either someone messed with this code or all logic in the world has been lost\n Let''s hope its the first one\n')
end

% Update class:
obj.DMpos = DMpos;
obj.ONorOFF = onoroff;

end % That's all Folks
