function plotLTspiceSchematic(fn, parser, positions)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
    
     %% ASC-like format
    %If an asc file is available, try to arrange the graph to match the
    %layout of the components in the schematic

    if isempty(parser.ascfn)
        schemFile = join(split(parser.sourcefn, ".net"));
        schemFile = [schemFile{:}(1:end-1) '.asc'];
    else 
        schemFile = parser.ascfn;
    end

    comps = parser.origComponents;
    compNames = {comps.Name};

    if exist(schemFile,"file") || (~exist('positions','var') || isempty(positions) )
        fID = fopen(schemFile);
%         schem = textscan(fID, '%s %d %d')
        schem = char(fread(fID)');
        fclose(fID);

        schem = splitlines(schem);
        schem(startsWith(schem,'WINDOW')) = [];
        wires = schem(startsWith(schem,'WIRE'));
        schem(startsWith(schem,'WIRE')) = [];
        componentStart = startsWith(schem,'SYMBOL');
        componentName = startsWith(schem,'SYMATTR InstName');

        locLines = schem(componentStart);
        nameLines = schem(componentName);

        names = strrep(nameLines, 'SYMATTR InstName ', '');
        locs = regexp(locLines, '(?<=SYMBOL [a-z]+ ).*', 'match');
        locs = split([locs{:}],' ');


        
        assocComponents = cell(numel(comps),1);
        xloc = zeros(numel(comps),1);
        yloc = zeros(numel(comps),1);
        rot = zeros(numel(comps),1);
        for i = 1:length(names)
            inds = endsWith(compNames, names{i});
            assocComponents(inds) = {names{i}};
            xloc(inds) = str2num(locs{1,i,1});
            yloc(inds) = str2num(locs{1,i,2});
            rot(inds) = str2num(locs{1,i,3}(2:end));
        end

        unidentified = cellfun(@isempty,assocComponents);

        positions = [xloc,yloc,rot];
    end

    %% Plot elements
    f = figure(fn); hold on; axis equal

    for i = 1:numel(comps)   
        if comps(i).Type == 'M'
            sf = 100/100;
            offset = [48,-50]; 
            compType = 'nmos';
        elseif comps(i).Type == 'L'
            sf = 2*82/100;
            offset = [16,-55];
            compType = 'ind';
        elseif comps(i).Type == 'C'
            sf = 70/100;
            offset = [16,-33];
            compType = 'cap';
        elseif comps(i).Type == 'R'
            sf = 80/100;
            offset = [16,-56]; 
            compType = 'res';
         elseif comps(i).Type == 'V'
            sf = 80/100;
            offset = [0,-56]; 
            compType = 'vsrc';
        end
        th = atan2(offset(:,2),offset(:,1)) -   positions(i,3)/180*pi;
        offset = sqrt(sum(offset.^2,2)).*[cos(th) sin(th)];
        plotCircuitElement(f,compType, [positions(i,1), -positions(i,2)]+offset, 100*sf,  positions(i,3), 1)
    end

    for i = 1:length(wires)
         coords = extract(wires{i},(""|"-") + digitsPattern);
         coords = cellfun(@str2double,coords);
         plot(coords([1 3]), -coords([2 4]),'k','LineWidth',1)
    end

    ax = f.Children(1);

    ax.XTick = [];
    ax.YTick = [];

        
        % %% Test
        % for i = find(unidentified)'
        %     ind = [];
        %     if startsWith(fullTable.Name{i},'V') && any(strcmp(names,['I' fullTable.Name{i}(2:end)]))
        %         ind = find(strcmp(names,['I' fullTable.Name{i}(2:end)]));
        %     elseif startsWith(fullTable.Name{i},'I') && any(strcmp(names,['V' fullTable.Name{i}(2:end)]))
        %         ind = find(strcmp(names,['V' fullTable.Name{i}(2:end)]));
        %     elseif length(split(fullTable.Name{i}, '_')) >= 3
        %         tryName = split(fullTable.Name{i}, '_');
        %         tryName = tryName{end-1};
        %         ind = find(strcmp(names,tryName));
        %     end
        % 
        %     if strcmpi(fullTable.Name{i},'W')
        %         xloc(i) = nan;
        %         yloc(i) = nan;
        %     else
        %         assocComponents(i) = {names{ind}};
        %         xloc(i) = str2num(locs{1,ind,1});
        %         yloc(i) = str2num(locs{1,ind,2});
        %         rot(i) = str2num(locs{1,ind,3}(2:end));
        %     end
        % end
        % 
        % fullTable = horzcat(fullTable, table(assocComponents, xloc, yloc, rot));

        % for i = 1:numel(assocComponents)
        %     if ~isempty(assocComponents(i))
        %         xData(i) = 
        % compIndices = ~strcmpi(fullTable.Name,'W');
        % wireIndices = strcmpi(fullTable.Name,'W');
        % compxData(compIndices) = fullTable.xloc(compIndices) + 5*cos(fullTable.rot(compIndices)/180*pi);
        % compyData(compIndices) = fullTable.yloc(compIndices) + 5*sin(fullTable.rot(compIndices)/180*pi);
        % 

        % uniqueNodes = unique(fullTable.EndNodes(:))';
        % for i = uniqueNodes
        %     connectedComponents = fullTable.EndNodes(:,1) == uniqueNodes(i) | fullTable.EndNodes(:,2) == uniqueNodes(i);
        %     if ~all(strcmpi(fullTable.Name(connectedComponents),'W'))
        %         compLoc = find(~strcmpi(fullTable.Name,'W') & connectedComponents,1);
        %         if fullTable.EndNodes(compLoc,1) == uniqueNodes(i)
        %             LR = 1;
        %         else
        %             LR = 0;
        %         end
        %         xData(i) = fullTable.xloc(compLoc) + LR*100*sin(-fullTable.rot(compLoc)/180*pi);
        %         yData(i) = fullTable.yloc(compLoc) + LR*100*cos(fullTable.rot(compLoc)/180*pi);
        % 
        %     else
        % 
        %         secondNodes = unique(fullTable.EndNodes(connectedComponents,:));
        %         [~,matchFirst] = intersect(fullTable.EndNodes(:,1) , secondNodes);
        %         [~,matchSecond] = intersect(fullTable.EndNodes(:,2) , secondNodes);
        %         secondConnections =  zeros(1,size(fullTable.EndNodes,1));
        %         secondConnections(matchFirst) = 1;
        %         secondConnections(matchSecond) = 1;
        %         connectedNonWires = secondConnections' & ~strcmpi(fullTable.Name,'W');
        %         xData(i) = mean(fullTable.xloc(connectedNonWires));
        %         yData(i) = mean(fullTable.yloc(connectedNonWires));
        %         if isnan(xData(i))
        %             xData(i) = 0;
        %             yData(i) = 0;
        %         end
        %     end
        % end
        % 
        % P.XData = xData;
        % P.YData = -yData;
    % 
    % end
end