function plotCircuitGraph2(obj)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    % if nargin == 1
    %     subgraph = obj.NL;
    % end
    % G = graph(subgraph(:,2)', subgraph(:,3)', subgraph(:,4)');
    % G.Edges.Name = {obj.components(G.Edges.Weight).Name}';
    % 
    % G2 = table();
    % nextPow10 = 10^ceil(log10(size(G.Edges,1)));
    % for i = size(G.Edges,1):-1:1
    %     segMentedComp = [G.Edges(i,:);
    %         G.Edges(i,:);
    %         G.Edges(i,:)];
    %     segMentedComp.Name(1) = {'W'};
    %     segMentedComp.Name(3) = {'W'};
    %     segMentedComp.EndNodes(1,2) =  nextPow10 + 2*i-1;
    %     segMentedComp.EndNodes(2,1) =  nextPow10 + 2*i-1;
    %     segMentedComp.EndNodes(2,2) = nextPow10 + 2*i;
    %     segMentedComp.EndNodes(3,1) =  nextPow10 + 2*i;
    %     G2 = [G2; segMentedComp];
    % end
    % 
    % curNodes = sort(unique(G2.EndNodes));
    % remap = dictionary(curNodes, [1:numel(curNodes)]');
    % G2.EndNodes = remap(G2.EndNodes);
    % 
    % G = graph(G2);
    % 
    % P = plot(G);
    % P.EdgeLabel = G.Edges.Name;


    %% ASC-like format
    %If an asc file is available, try to arrange the graph to match the
    %layout of the components in the schematic

    if isempty(obj.ascfn)
        schemFile = join(split(obj.sourcefn, ".net"));
        schemFile = [schemFile{:}(1:end-1) '.asc'];
    else 
        schemFile = obj.ascfn;
    end

    if exist(schemFile,"file")
        fID = fopen(schemFile);
%         schem = textscan(fID, '%s %d %d')
        schem = char(fread(fID)');
        fclose(fID);

        schem = splitlines(schem);
        schem(startsWith(schem,'WINDOW')) = [];
        schem(startsWith(schem,'WIRE')) = [];
        componentStart = startsWith(schem,'SYMBOL');

        locLines = schem(componentStart);
        nameLines = schem([0==1; componentStart(1:end-1)]);

        names = strrep(nameLines, 'SYMATTR InstName ', '');
        locs = regexp(locLines, '(?<=SYMBOL [a-z]+ ).*', 'match');
        locs = split([locs{:}],' ');


        fullTable = G.Edges;
        assocComponents = cell(length(fullTable.Name),1);
        xloc = zeros(length(fullTable.Name),1);
        yloc = zeros(length(fullTable.Name),1);
        rot = zeros(length(fullTable.Name),1);
        for i = 1:length(names)
            inds = endsWith(fullTable.Name, names{i});
            assocComponents(inds) = {names{i}};
            xloc(inds) = str2num(locs{1,i,1});
            yloc(inds) = str2num(locs{1,i,2});
            rot(inds) = str2num(locs{1,i,3}(2:end));
        end

        unidentified = cellfun(@isempty,assocComponents);
        %% Test
        for i = find(unidentified)'
            ind = [];
            if startsWith(fullTable.Name{i},'V') && any(strcmp(names,['I' fullTable.Name{i}(2:end)]))
                ind = find(strcmp(names,['I' fullTable.Name{i}(2:end)]));
            elseif startsWith(fullTable.Name{i},'I') && any(strcmp(names,['V' fullTable.Name{i}(2:end)]))
                ind = find(strcmp(names,['V' fullTable.Name{i}(2:end)]));
            elseif length(split(fullTable.Name{i}, '_')) >= 3
                tryName = split(fullTable.Name{i}, '_');
                tryName = tryName{end-1};
                ind = find(strcmp(names,tryName));
            end
            
            if strcmpi(fullTable.Name{i},'W')
                xloc(i) = nan;
                yloc(i) = nan;
            else
                assocComponents(i) = {names{ind}};
                xloc(i) = str2num(locs{1,ind,1});
                yloc(i) = str2num(locs{1,ind,2});
                rot(i) = str2num(locs{1,ind,3}(2:end));
            end
        end

        fullTable = horzcat(fullTable, table(assocComponents, xloc, yloc, rot));

        % for i = 1:numel(assocComponents)
        %     if ~isempty(assocComponents(i))
        %         xData(i) = 
        % compIndices = ~strcmpi(fullTable.Name,'W');
        % wireIndices = strcmpi(fullTable.Name,'W');
        % compxData(compIndices) = fullTable.xloc(compIndices) + 5*cos(fullTable.rot(compIndices)/180*pi);
        % compyData(compIndices) = fullTable.yloc(compIndices) + 5*sin(fullTable.rot(compIndices)/180*pi);
        % 

        uniqueNodes = unique(fullTable.EndNodes(:))';
        for i = uniqueNodes
            connectedComponents = fullTable.EndNodes(:,1) == uniqueNodes(i) | fullTable.EndNodes(:,2) == uniqueNodes(i);
            if ~all(strcmpi(fullTable.Name(connectedComponents),'W'))
                compLoc = find(~strcmpi(fullTable.Name,'W') & connectedComponents,1);
                if fullTable.EndNodes(compLoc,1) == uniqueNodes(i)
                    LR = 1;
                else
                    LR = 0;
                end
                xData(i) = fullTable.xloc(compLoc) + LR*100*sin(-fullTable.rot(compLoc)/180*pi);
                yData(i) = fullTable.yloc(compLoc) + LR*100*cos(fullTable.rot(compLoc)/180*pi);

            else
            
                secondNodes = unique(fullTable.EndNodes(connectedComponents,:));
                [~,matchFirst] = intersect(fullTable.EndNodes(:,1) , secondNodes);
                [~,matchSecond] = intersect(fullTable.EndNodes(:,2) , secondNodes);
                secondConnections =  zeros(1,size(fullTable.EndNodes,1));
                secondConnections(matchFirst) = 1;
                secondConnections(matchSecond) = 1;
                connectedNonWires = secondConnections' & ~strcmpi(fullTable.Name,'W');
                xData(i) = mean(fullTable.xloc(connectedNonWires));
                yData(i) = mean(fullTable.yloc(connectedNonWires));
                if isnan(xData(i))
                    xData(i) = 0;
                    yData(i) = 0;
                end
            end
        end

        P.XData = xData;
        P.YData = -yData;

    end



end