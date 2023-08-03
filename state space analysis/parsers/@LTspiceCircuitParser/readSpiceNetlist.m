function readSpiceNetlist(obj,filename)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    fid=fopen(filename);
    text = char(fread(fid)'); 
    fclose(fid);

    floc = dir(which(filename));
    if isempty(floc)
        obj.sourcefdate = '0';
    else
        obj.sourcefdate = floc.date;
    end
    
    lines = splitlines(text);

    ind4 = @(x) {x(6:end)};
    ind6 = @(x) {x(8:end)};

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
    for i = 1:length(componentsLines)
        newComponent = obj.parseSpiceComponent(componentsLines{i});
        for j = 1:length(newComponent) % second loop in case parsing dembedded parameters
            if strcmp(newComponent(j).Type, 'V') || strcmp(newComponent(j).Type, 'I')
                sources = [sources, newComponent(j)];
            elseif strcmp(newComponent(j).Type, 'M') || strcmp(newComponent(j).Type, 'D')
                switches = [switches, newComponent(j)];
            elseif strcmp(newComponent(j).Type, 'K') 
                coupling = [coupling, newComponent(j)];
            else
                passives = [passives, newComponent(j)];
            end
        end
    end

    components = [passives, sources, switches];


    obj.origComponents = components;
    obj.netListDirectives = coupling;

    obj.topology.switchLabels = {switches.Name};


end