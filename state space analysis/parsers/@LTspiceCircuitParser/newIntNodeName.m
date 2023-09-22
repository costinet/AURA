function [nodeName] = newIntNodeName(obj,component,node, num)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    if nargin == 3
        strnum = '';
    else
        strnum = ['_' num2str(num)];
    end

    nodeName = [component.Nodes{node} 'm_' component.Name strnum];

    %Make sure there won't be an accidental collision with a node name from
    %the original file
    i=2;
    if ~isempty(obj.origComponents)
        while any(strcmp([obj.origComponents.Nodes], nodeName))
            if ~exist('baseNodeName','var')
                baseNodeName = nodeName;
            end
            nodeName = [baseNodeName '_' num2str(i)];
            i = i+1;
        end
    end
end