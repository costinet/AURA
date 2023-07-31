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
    for i = 1:length(componentsLines)
        newComponent = obj.parseSpiceComponent(componentsLines{i});
        if strcmp(newComponent.Type, 'V') || strcmp(newComponent.Type, 'I')
            sources = [sources, newComponent];
        elseif strcmp(newComponent.Type, 'M') || strcmp(newComponent.Type, 'D')
            switches = [switches, newComponent];
        else
            passives = [passives, newComponent];
        end
    end

    components = [passives, sources, switches];

    %% Interpret
    % Nodes
    nodes = unique([components(:).Nodes]);
    nodes = circshift(nodes, -1*(find(strcmp(nodes,'0'))-1)); % make '0' the first node
    nodeMap = dictionary(nodes, 0:length(nodes)-1);

    % Component Types
    typeMap = dictionary({'V','BV','MV','C','R','L','MI','BI','I'}, 1:9);

    %% for later? Numerical nodes corresponding to components
    numNodes = reshape(nodeMap([components.Nodes]), [2,length(components),])';
%     typeOrder = typeMap({components.Type});  % for his formatting, this
%     has to be done after replacing switches with Rs and Cs
                                                

    %% Compatibility with Jared's formatting
    FETs = [components.Type] == 'M';
    diodes = [components.Type] == 'D';
    switches = FETs | diodes;
    numSwitches = sum(switches);
    numDiodes = sum(diodes);

    obj.Switch_Names = {components(switches).Name}';
    obj.Switch_Resistors = strcat(obj.Switch_Names, '_R');
%     obj.Switch_Resistor_Values = cell2mat({components(switches).paramVals(1:3)}');
%     obj.Switch_Resistor_Values = [obj.Switch_Resistor_Values(:,2) ...
%         obj.Switch_Resistor_Values(:,1) obj.Switch_Resistor_Values(:,2)];
    obj.Switch_Resistor_Values = [];

    components = [components, components(switches)];
%     [components(end-numSwitches+1:end).Type] = deal('C');
    locs = find(switches);
    for i = 1:numSwitches 
        loc = locs(i);

        obj.Switch_Resistor_Values(i,:) = [
            components(loc).paramVals(strcmp(components(loc).paramNames, 'Roff')), ...
            components(loc).paramVals(strcmp(components(loc).paramNames, 'Rds') | strcmp(components(loc).paramNames, 'Ron')), ...
            components(loc).paramVals(strcmp(components(loc).paramNames, 'Roff'))
            ];


        components(loc).Name =  [strcat(components(loc).Name,'_C')];
        components(loc).Type = 'C';
        components(loc).paramNames = {'C'};
        components(loc).paramVals = components(loc).paramVals(3);

        
        components(end-i+1).Name = [strcat(components(end-i+1).Name,'_R')];
        components(end-i+1).Type = 'R';
        components(end-i+1).paramNames = {'R'};
        components(end-i+1).paramVals = components(end-i+1).paramVals(2);
    end

    obj.components = components;

    [~,IA,~] = intersect({components.Name}, {sources.Name});
    obj.Element_Properties = [{components.Name}' {components.paramVals}'];
    obj.Element_Properties(IA,:) = [];


end