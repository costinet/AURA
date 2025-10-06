function [components] = findShorts(obj, components)
%findShorts removed any shorted components (R->0 or L->0) and updates nodes
%on other elements to remove redundant node labels
%    Predominately needed when updating but not reparsing.  origComponents
%    will retain the element.

    %Need to neglect depedent sources because they have multiple items in
    %paramVals
    compVals = {components.paramVals};
    compVals(strcmpi({components.Type},'E') | strcmpi({components.Type},'F') | ...
        strcmpi({components.Type},'G') | strcmpi({components.Type},'H')) = {0};

    shorts = (strcmpi({components.Type},'L') & cell2mat(compVals) == 0)...
        | ... % or
        (strcmpi({components.Type},'R') & cell2mat(compVals) == 0);
    
    for i = find(shorts)
        shortedNodes = components(i).Nodes;
        for j = find(~shorts)
            % replace redundant node in every other component
            components(j).Nodes(strcmp(components(j).Nodes,shortedNodes(1))) = shortedNodes(2);
        end
    end
    components(shorts) = [];


end