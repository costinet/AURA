function [ALL,NL,NL3,incidence] = read_file(filename)
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
% Capacitor: 2
% Resistor: 3
% Inductor: 4
% Current Source: 5
% Diode: 6
% FET: 7
% Dependent Voltage Source: 8
% Dependent Current Source: 9
% 
% 2nd and 3rd column provide the nodes the corressponding element described
% in 1st column
% 
% Incidence: Provides and incidence matrix of the elements in the circuit
%
%
% filename is a character vector or string scalar that is the name of the
% NETlist to be opened
%
% Example: 'BoostConverter.net'
% 





if exist(filename,'file') == 2
    
    fid=fopen(filename);
    tline = [];
    i=0;
    
    % Open file and read into cell array:
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
    error(['Error: ',filename,' Does Not Exist'])
end

% Key for elements:
% Value     1   2   3   4   5   6   7   8
pattern = {'V','C','R','L','I','D','M','B'};

numV = 1;
% BV = 2;
numC = 2;
numR = 3;
numL = 4;
numI = 5;
% BI = 6;
numD = 6;
numM = 7;
numB = 8;


m = 0;

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
            for j = 1:1:length(Net) 
                ALL(i,j) = Net(j);  % Add contents of row to all array
                NL(1+m,j) = Net(j); % Add contents of row to character netlist only array array
                NL2(1+m,:) = [k;Net(2);Net(3);1+m]; % Change letter to number add contents to numerical netlist array
                if k == 7 % Fix FET Nodes
                    NL2(1+m,3)= Net(5); % Re assign nodes for fets to have drain and source in rows 2 and 3
                end
                
                %%%%%%%%% Start New
                
                if k == 8 && j==4 
                    
                    temp = Net{j};
                    
                    if strcmp(temp(1),'I')
                        NL2(1+m,1)={9};
                    end
                    
                end
                
                %%%%%%% End new
                
            end
            m = m+1;
        end
    end
end
    
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

% Initally sort rows into normal tree priorety 
% NL3 = sortrows(NL3,1);

incidence = zeros(max(new_list),length(NL3)); % initialize incidence matrix
% Find Incidence matrix
for i = 1:1:length(NL3)
    incidence(NL3(i,2),i) = 1;
    incidence(NL3(i,3),i) = -1;
end
    
if isempty(NL3) || isempty(incidence)
    error('WARNING: No valid elements found in circuit');
end


