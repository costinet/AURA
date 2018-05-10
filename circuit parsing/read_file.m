function [ALL,NL,NL3,K] = read_file(filename)
% 
% The code reads the NETlist file from LTSpice to be used in finding ABCD
%
% [NLraw,ALL,NL,NL3,incidence] = read_file(filename) returns:
% 
% ALL: Gives the entire content of the NETlist parsed by lines and spaces
% in a cell array
% 
% NL: Gives only the components of the circuit parsed by lines and spaces
% in a cell arrary
% 
% NL3: 1st column provides a list of the components in the circuit
% corresponding the the following legend
% Voltage Source: 1
% Dependent Voltage Source (Transformer): 2
% Measurement Voltage Source: 3
% Capacitor: 4
% Resistor: 5
% Inductor: 6
% Measurement Current Source: 7
% Dependent Current Source (Transformer): 8
% Current Source: 9
% Diode: 10
% FET: 11
% 
% 2nd and 3rd column provide the nodes the corressponding element described
% in 1st column
%
%
% filename is a character vector or string scalar that is the name of the
% NETlist to be opened
%
% Example: 'BoostConverter.net'
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

 
%% Open File

if exist(filename,'file') == 2 % if the file exists
    
    fid=fopen(filename); % open file
    tline = [];
    i=0;
    
    % Read file into cell array:
    while ~feof(fid)
        tline = fgetl(fid);
        string = cellstr(tline);
        NLraw(1+i)=string;
        i=i+1;
    end
    
    % Check to ensure the file closed:
    check_close = fclose(fid);
    
    % If it didn't close display exit and display warning
    if check_close == -1
        error(['Error: ',filename,' Not Properly Closed'])
    end
    
else
    % If file does not exist display error
    error(['Error: ',filename,' Does Not Exist'])
end

% Key for elements:
% Value     1   2   3   4   5   6   7   8
pattern = {'V','C','R','L','I','D','M','K'};

numV = 1;
numC = 2;
numR = 3;
numL = 4;
numI = 5;
numD = 6;
numM = 7;
numK = 8;


m = 0;
IndtoTrans = [];
MutInd = [];


%% Parse data from file

for i = 1:1:length(NLraw) % Step though netlist
    cell = NLraw(i); % Take each row of the netlist
    chars = char(cell); % Convert to chars
    Net = strsplit(chars); % Split up chars by spaces into cells
    for k = 1:1:length(pattern) % Step through pattern of components
        emptyish = strfind(Net{1},pattern{k}); % find if component was in cell
        if isempty(emptyish) % if there was no componenets
            % Add contents of row to all array 
            for j = 1:1:length(Net)
                ALL(i,j) = Net(j);
            end
        else
            if strcmp(pattern{k},'K') % Check to see if there is a transformer
                IndtoTrans = Net(2:end-1); % If true then get inductors that are a part of transformer
                MutInd = Net(end); % Get mutual inductance value for transformer
                m = m - 1;
            else
                for j = 1:1:length(Net)
                    ALL(i,j) = Net(j);  % Add contents of row to all array
                    NL(1+m,j) = Net(j); % Add contents of row to character netlist only array array
                    NL2(1+m,:) = [k;Net(2);Net(3);1+m]; % Change letter to number add contents to numerical netlist array
                    if k == numM % Fix FET Nodes
                        NL2(1+m,3)= Net(5); % Re assign nodes for fets to have drain and source in rows 2 and 3
                    end
                    
                end
            end
            m = m+1;
        end
    end
end

%% Set 0 node to be ground

List = NL2(:,2:3); % list contains only the nodes
gndpt = cellfun(@(s) contains('0', s), List); % find points where there is a ground
gndpt = [gndpt(:,1);gndpt(:,2)];    % Put all values in one column
[~,~,new_list] = unique(List,'stable'); % Assign node numbers
lil=find((gndpt.*new_list)~=0,1); % Find index assigned to 0
replaced=new_list(lil); % Find value assigned to zero
for i = 1:1:length(new_list)
    if new_list(i)==1
        new_list(i) = replaced;
    end
    if 0 ~= gndpt(i) * new_list(i)% arrary of
        new_list(i) = 1;
    end
end
num = length(new_list)/2;

% Fix output in the form [component type, node1, node2]
NL3 = [cell2mat(NL2(:,1)),new_list(1:num,:),new_list(num+1:length(new_list),:),cell2mat(NL2(:,4))]; 

%{
% Initally sort rows into normal tree priorety 
% NL3 = sortrows(NL3,1);

incidence = zeros(max(new_list),length(NL3)); % initialize incidence matrix
% Find Incidence matrix
for i = 1:1:length(NL3)
    incidence(NL3(i,2),i) = 1;
    incidence(NL3(i,3),i) = -1;
end
  
  %}
if isempty(NL3) %|| isempty(incidence)
    error('WARNING: No valid elements found in circuit');
end

K = [];



%% Transformers and Measurements

% Set new index
% Name:       V BV MV  C  R  L MI BI  I  D  M
% Identifier: 1  2  3  4  5  6  7  8  9 10 11

% Set New Index code
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

% Implements index above:
NL3(:,1)=NL3(:,1)+(NL3(:,1)>4);
NL3(:,1)=NL3(:,1)+(NL3(:,1)>4);
NL3(:,1)=NL3(:,1)+(NL3(:,1)>1);
NL3(:,1)=NL3(:,1)+(NL3(:,1)>1);


% Find location of Transfomer elements in NL
D = size(IndtoTrans);
position = zeros(D);
for i = 1:1:D(1)
    for j = 1:1:D(2)
        found=strcmp(NL(:,1),IndtoTrans(i,j));
        position(i,j)=find(found,1);
    end
end
        
[NL3,NL] = addmeasure(NL3,NL);
% Find and set transformers
if ~isempty(IndtoTrans) || ~isempty(MutInd)
[NL,NL3,K]=transformers(NL,NL3,IndtoTrans,MutInd,position);
end


J = 89704895734895;

