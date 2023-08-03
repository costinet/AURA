function [A,B,C,D,I] = solveStateSpaceRepresentation(obj)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     read_file_num(obj) % needed to set obj.NL
%     [switches]=obj.findDM;
%     obj.ON_States = cell(length(switches),1);
%     obj.OFF_States = cell(length(switches),1);
%     [obj.NewNL,obj.NewNLnets,~]=obj.Single_states_D(1,1,switches);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Needed anymore?
    sources = obj.components(strcmp({obj.components.Type},'V') | strcmp({obj.components.Type},'I'));
    [~,IA,~] = intersect({obj.components.Name}, {sources.Name});
    obj.Element_Properties = [{obj.components.Name}' {obj.components.paramVals}'];
    obj.Element_Properties(IA,:) = [];
    obj.Component_Values = obj.Element_Properties;


    %% Interpret
    % Give nodes numeric identifiers
    nodes = unique([obj.components(:).Nodes]);
    nodes = circshift(nodes, -1*(find(strcmp(nodes,'0'))-1)); % make '0' the first node
    nodeMap = dictionary(nodes, 1:length(nodes));
    numNodes = reshape(nodeMap([obj.components.Nodes]), [2,length(obj.components),])';

    %% for later? Numerical nodes corresponding to components
    typeMap = dictionary({'V','BV', 'MV','C','R','L','MI','BI','I', 'E', 'F'}, [1:9, 2, 8]);
    typeOrder = typeMap({obj.components.Type})';  

    obj.NL = [typeOrder, numNodes , transpose(1:length(typeOrder))];
    obj.NLnets = [{obj.components.Name}' {obj.components.Nodes}' ];

    %% Graph backend?
%     G = graph(numNodes(:,1)', numNodes(:,2)', 1:size(numNodes,1));
%     G.Edges.Name = {components(G.Edges.Weight).Name}';
% 
%     P = plot(G);
%     P.EdgeLabel = G.Edges.Name;



    obj.cutset_loop_num();
    [A,B,C,D,I] = obj.nodeloop_num(obj.NL,obj.NLnets);
end