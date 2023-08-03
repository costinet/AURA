function linearizeCircuitModel(obj)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

    components = obj.origComponents;
    
    %% Replace FETs with Linear equivalent circuit
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

    sources = components(strcmp({components.Type},'V') | strcmp({components.Type},'I'));
    [~,IA,~] = intersect({components.Name}, {sources.Name});
    obj.Element_Properties = [{components.Name}' {components.paramVals}'];
    obj.Element_Properties(IA,:) = [];



% % %     %% Interpret
% % %     % Nodes
% % %     nodes = unique([components(:).Nodes]);
% % %     nodes = circshift(nodes, -1*(find(strcmp(nodes,'0'))-1)); % make '0' the first node
% % %     nodeMap = dictionary(nodes, 1:length(nodes));
% % % 
% % %     % Component Types
% % %     typeMap = dictionary({'V','BV','MV','C','R','L','MI','BI','I'}, 1:9);
% % % 
% % %     %% for later? Numerical nodes corresponding to components
% % %     numNodes = reshape(nodeMap([components.Nodes]), [2,length(components),])';
% % % %     typeOrder = typeMap({components.Type});  % for his formatting, this
% % % %     has to be done after replacing switches with Rs and Cs
% % % 
% % %     %% Graph backend?
% % % %     G = graph(numNodes(:,1)', numNodes(:,2)', 1:size(numNodes,1));
% % % %     G.Edges.Name = {components(G.Edges.Weight).Name}';
% % % % 
% % % %     P = plot(G);
% % % %     P.EdgeLabel = G.Edges.Name;
end