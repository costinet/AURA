function plotCircuitGraph(obj, subgraph)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    if nargin == 1
        subgraph = obj.NL;
    end
    G = graph(subgraph(:,2)', subgraph(:,3)', subgraph(:,4)');
    G.Edges.Name = {obj.components(G.Edges.Weight).Name}';

    P = plot(G);
    P.EdgeLabel = G.Edges.Name;


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

            assocComponents(i) = {names{ind}};
            xloc(i) = str2num(locs{1,ind,1});
            yloc(i) = str2num(locs{1,ind,2});
            rot(i) = str2num(locs{1,ind,3}(2:end));
        end

        fullTable = horzcat(fullTable, table(assocComponents, xloc, yloc, rot));

        for i = unique(fullTable.EndNodes(:))'
            connectedComponents = fullTable.EndNodes(:,1) == i | fullTable.EndNodes(:,2) == i;
            xData(i) = mean(fullTable.xloc(connectedComponents));
            yData(i) = mean(fullTable.yloc(connectedComponents));
        end

        P.XData = xData;
        P.YData = -yData;

    end



end