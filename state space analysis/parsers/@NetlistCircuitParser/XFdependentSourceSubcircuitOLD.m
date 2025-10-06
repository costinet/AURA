function components = XFdependentSourceSubcircuitOLD(obj, directive, components)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

    
    for i = 1:length(directive.paramNames)
        coupledIndLocs(i) = find(strcmp({components.Name}, directive.paramNames{i}));
        coupledInductors(i) = components(coupledIndLocs(i));
    end
    try
        k = eval(directive.paramVals);
    catch
         paramVal = directive.paramVals;
         paramVal = strrep(paramVal,'{','');
         paramVal = strrep(paramVal,'}','');
         k = evalin('base',paramVal);
    end

    if k>=1
        %% Use ideal transformer model
        % model has one current source and all other ports will be voltage
        % sources.  The following block tries to find a suitable place to
        % put the current source where it won't results in L-I loops

        %% Attempt to find a low-impedance place to put the current source
        nodes = unique([components(:).Nodes]);
        nodes = circshift(nodes, -1*(find(strcmp(nodes,'0'))-1)); % make '0' the first node
        nodeMap = dictionary(nodes, 1:length(nodes));
        numNodes = reshape(nodeMap([components.Nodes]), [2,length(components),])';

        possibleLocs = zeros(length(coupledIndLocs),1);

        % This only looks at the first-order connections.  Loops would be
        % better.
        for i = 1:length(coupledIndLocs)
            for j = 1:2
                Node = numNodes(coupledIndLocs(i),j);
                connectedComponents = [find(numNodes(:,1) == Node);  find(numNodes(:,2) == Node)];
                connectedComponents = unique(connectedComponents);
                connectedComponents(connectedComponents == coupledIndLocs(i)) = [];

                connectedComponentsTypes = {components(connectedComponents).Type};
                if isempty(connectedComponentsTypes)
                    possibleLocs(i) = -2;
                elseif all(strcmp(connectedComponentsTypes, "L") | strcmp(connectedComponentsTypes, "I")  | ...
                        strcmp(connectedComponentsTypes, "Im")) 
                    possibleLocs(i) = -1;
                end
            end
        end

        Vi = find(possibleLocs == max(possibleLocs),1);

        % Vi = find(possibleLocs>=0,1,'first');
        % 
        % %% If we couldn't find a good location, just guess
        % if(isempty(Vi))
        %     Vi = 1; % for the moment, assume the first inductor is the voltage reference.  
        %             % later, should look at impedance seen by each coil
        % end

        %% Replace sources
        N(Vi) = 1;
        Nexpr{Vi} = '1';
        LV = coupledInductors(Vi).paramVals(strcmp(coupledInductors(Vi).paramNames,'L'));
        isource = {};
        isource.Name = ['F_' coupledInductors(Vi).Name];
        isource.Type = 'F';
        isource.Nodes = coupledInductors(Vi).Nodes;
        isource.paramNames = {};
        isource.paramVals = {};
        isource.paramExpressions = {};

        Lmin = 0;

        for i = setdiff(1:length(coupledInductors),Vi)
            Li = coupledInductors(i).paramVals(strcmp(coupledInductors(i).paramNames,'L'));
            N(i) = sqrt(Li/LV);
            Nexpr{i} = ['sqrt(' coupledInductors(i).paramExpressions(strcmp(coupledInductors(i).paramNames,'L')) ...
                '/' coupledInductors(Vi).paramExpressions(strcmp(coupledInductors(Vi).paramNames,'L')) ')'];
            
            Vsource = {};
            Vsource.Name = ['E_' coupledInductors(i).Name];
            Vsource.Type = 'E';
            Vsource.Nodes = [coupledInductors(i).Nodes];
            Vsource.paramNames = {'Inam', 'gain'};   %This is not valid LTspice, but it is easier to make it a two-node component
            Vsource.paramVals = {isource.Name, N(i)};
            Vsource.paramExpressions = {isource.Name, Nexpr{i} };

            isource.paramNames = [isource.paramNames, ...
                {['Vnam' num2str(i)], ['gain' num2str(i)]}];
            isource.paramVals = [isource.paramVals, ...
                {Vsource.Name, -N(i)}];
            isource.paramExpressions = {isource.paramExpressions; {isource.Name, ['-' Nexpr{i}]} };

            if Lmin == 0
                Lmin = 1;
                %Leave one Lm in (save for later to not ruin indexing
                Lm = components(coupledIndLocs(i));  
            end
            components(coupledIndLocs(i)) = Vsource;    %Replace all secondaries

        end


        components(coupledIndLocs(Vi)) = isource;
        if ~(k>1)
            components(end+1) = Lm;
            % components = [components(1:NVi-1) isource components(NVi+1:end)]; 
        end
    else
        error('coupled inductors with k < 1 not currently supported');
    end

end