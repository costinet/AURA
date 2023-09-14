function [nodeName] = newIntNodeName(obj,component,node, num)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    if nargin == 3
        strnum = '';
    else
        strnum = ['_' num2str(num)];
    end

    nodeName = [component.Nodes{node} 'm_' component.Name strnum];
end