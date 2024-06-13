function readLTspiceSchematic(obj)
%readLTspiceSchematic read *.asc file to 
%   Detailed explanation goes here

    %% ASC-like format
    %If an asc file is available, try to arrange the graph to match the
    %layout of the components in the schematic

    if isempty(obj.ascfn)
        schemFile = join(split(obj.sourcefn, ".net"));
        schemFile = [schemFile{:}(1:end-1) '.asc'];
    else 
        schemFile = obj.ascfn;
    end

    comps = obj.origComponents;
    compNames = {comps.Name};

    % Replace measurement sources back to V/I so they can be found
    compNames(strcmp({comps.Type},'Vm')) = ...
        cellfun(@(x)['V' x(2:end)], compNames(strcmp({comps.Type},'Vm')), 'UniformOutput', false);
    compNames(strcmp({comps.Type},'Im')) = ...
        cellfun(@(x)['I' x(2:end)], compNames(strcmp({comps.Type},'Im')), 'UniformOutput', false);

    if exist(schemFile,"file")
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
            if any(inds)
                assocComponents(inds) = {names{i}};
                xloc(inds) = str2num(locs{1,i,1});
                yloc(inds) = str2num(locs{1,i,2});
                rot(inds) = str2num(locs{1,i,3}(2:end));
            end
        end

        unidentified = cellfun(@isempty,assocComponents);

        positions = [xloc,yloc,rot];
    else
        return;
    end

    schemPos = zeros(numel(comps),4); %[X, Y, rot, SF]
    for i = 1:numel(comps)   
        if comps(i).Type == 'M'
            sf = 100/100;
            offset = [48,-50]; 
            % compType = 'nmos';
        elseif comps(i).Type == 'L'
            sf = 82/100;
            offset = [16,-55];
            % compType = 'ind';
        elseif comps(i).Type == 'C'
            sf = 70/100;
            offset = [16,-33];
            % compType = 'cap';
        elseif comps(i).Type == 'R'
            sf = 80/100;
            offset = [16,-56]; 
            % compType = 'res';
         elseif startsWith(comps(i).Type,'V')
            sf = 80/100;
            offset = [0,-56]; 
            % compType = 'vsrc';
        end
        th = atan2(offset(:,2),offset(:,1)) -   positions(i,3)/180*pi;
        offset = sqrt(sum(offset.^2,2)).*[cos(th) sin(th)];
        schemPos(i,:) = [[positions(i,1), -positions(i,2)]+offset, positions(i,3), 100*sf];
    end

    obj.schemPositions = schemPos;

    wireCoords = zeros(length(wires),4);
    for i = 1:length(wires)
         coords = extract(wires{i},(""|"-") + digitsPattern);
         wireCoords(i,:) = cellfun(@str2double,coords); 
    end

    obj.schemWires = wireCoords;

end