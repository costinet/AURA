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

end