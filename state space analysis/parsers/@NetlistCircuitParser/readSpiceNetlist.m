function readSpiceNetlist(obj,filename)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    if exist(filename, 'file') == 2
        fid=fopen(filename);
        text = char(fread(fid)'); 
        fclose(fid);
        obj.sourceType = 'file';
        obj.sourcefn = filename;

        floc = dir(which(filename));
        if isempty(floc)
            obj.sourcefdate = '0';
        else
            obj.sourcefdate = floc.date;
        end
    else
        text = filename;
        obj.sourceType = 'string';
    end
    
    lines = splitlines(text);

    ind4 = @(x) {x(6:end)};
    ind6 = @(x) {x(8:end)};

    obj.undefinedExpressions = {};

    libs = [cellfun(ind4, lines(startsWith(lines,'.lib')))];
    params = [cellfun(ind6, lines(startsWith(lines,'.param')))];
    models = [cellfun(ind6, lines(startsWith(lines,'.model')))];
    componentsLines = lines(~startsWith(lines,'.') & ~startsWith(lines,'*') & ~cellfun(@isempty,lines));
    
    %% Parameters
    % Handle params first so if anything is needed later, it will be available
    % from the base workstation
    for i = 1:length(params)
        obj.evalSpiceParams([params{i}])
    end
    
    %% Models and Libraries
    % Some of these are defaults, so we'll save them away and see if we
    % need them later rather than reading them all in
    obj.netlistLibraries = libs;
    obj.netlistModels = models;
    
    %% Components
    passives = [];
    sources = [];
    switches = [];
    coupling = [];
    meters = [];
    for i = 1:length(componentsLines)
        newComponent = obj.parseSpiceComponent(componentsLines{i});
        for j = 1:length(newComponent) % second loop in case parsing dembedded parameters
            if strcmp(newComponent(j).Type, 'V') || strcmp(newComponent(j).Type, 'I')
                sources = [sources, newComponent(j)];
            elseif strcmp(newComponent(j).Type, 'M') || strcmp(newComponent(j).Type, 'D')
                switches = [switches, newComponent(j)];
            elseif strcmp(newComponent(j).Type, 'K') 
                coupling = [coupling, newComponent(j)];
            elseif strcmp(newComponent(j).Type, 'Vm') || strcmp(newComponent(j).Type, 'Im')
                meters = [meters, newComponent(j)];
            else
                passives = [passives, newComponent(j)];
            end
        end
    end

    %% 
    components = [passives, sources, switches, meters];



    %% Check for zero-valued passives
    shorts = (strcmp({passives.Type},'L') & ([passives.paramVals] == 0) ) | ...
        (strcmp({passives.Type},'R') & ([passives.paramVals] == 0) );
    opens = (strcmp({passives.Type},'C') & ([passives.paramVals] == 0) ) | ...
        (strcmp({passives.Type},'L') & isinf([passives.paramVals]) ) | ...
        (strcmp({passives.Type},'R') & isinf([passives.paramVals]) );

    passives(opens) = [];

    sameNodes = {passives(shorts).Nodes}';

    passives(shorts) = [];

    components = [passives, sources, switches, meters];

    for i = 1:length(sameNodes)
        for j = 1:length(components)
            if any(strcmp(components(j).Nodes, sameNodes{i}{1}))
                components(j).Nodes{strcmp(components(j).Nodes, sameNodes{i}{1})} = ...
                    sameNodes{i}{2};
            end
        end
    end

    obj.origComponents = components;
    obj.netListDirectives = coupling;

    obj.topology.switchLabels = {switches.Name};

    try
        obj.readLTspiceSchematic;
    catch
        % Do nothing.  This function is only to gelp with graphical
        % representation
    end


end