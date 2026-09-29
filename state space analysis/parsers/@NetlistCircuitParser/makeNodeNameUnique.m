function [nodeNameO] = makeNodeNameUnique(obj,nodeName)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    obj circuitParser
    nodeName {mustBeText}
end

arguments (Output)
    nodeNameO {mustBeText}
end

nodeNameO = nodeName;
origNodes = [obj.origComponents.Nodes];
if ~isempty(obj.origComponents)
    origNodes = [obj.origComponents.Nodes];
else
    origNodes = {};
end

if ~isempty(obj.components)
    compNodes = [obj.components.Nodes];
else
    compNodes = {};
end

i=2;
while any(strcmp([origNodes compNodes], nodeName))
    if ~exist('baseNodeName','var')
        baseNodeName = nodeName;
    end
    nodeName = [baseNodeName '_' num2str(i)];
    i = i+1;
end
