function plotSchematic(obj,fn)
%plotSchematic generate visual schematic of original circuit
%   Detailed explanation goes here

    if nargin < 2
        f = fig;
    else
        f = figure(fn);
    end

    %% Plot elements
     hold on; axis equal

    comps = obj.origComponents;
    compNames = {comps.Name};

    % compNameDict = dictionary('M','nmos','L','ind','C','cap','R','res','V','vsrc','I','isrc','D','dio');

    if ~isempty(obj.schemPositions)
        for i = 1:numel(comps) 
            if ~any(strcmp({'Rser'},comps(i).paramNames))
                plotCircuitElement(f,comps(i).Type, obj.schemPositions(i,1:2),...
                    obj.schemPositions(i,4),  obj.schemPositions(i,3), 1);
                
                
                [Xoff,Yoff] = pol2cart(obj.schemPositions(i,3)*pi/180,10);
                % text(obj.schemPositions(i,1)+Xoff,obj.schemPositions(i,2)+Yoff, [compNames{i} newline comps(i).paramExpressions])
                
                if iscell(comps(i).paramExpressions) && numel(comps(i).paramExpressions) > 1
                    fulltext = strjoin({compNames{i}, ['\color[rgb]{.5 .5 .5} {\fontsize{8}' strjoin(comps(i).paramExpressions( ~cellfun(@isempty,comps(i).paramExpressions)), '; ' ) '}']},newline);
                else
                   fulltext = strjoin({compNames{i}, ['\color[rgb]{.5 .5 .5}  {\fontsize{8}' comps(i).paramExpressions '}']},newline);
                end
                % fulltext = strjoin({compNames{i}, strjoin(comps(i).paramExpressions( ~cellfun(@isempty,comps(i).paramExpressions)), '; ' )},newline);
                text(obj.schemPositions(i,1)+Xoff,obj.schemPositions(i,2)+Yoff, fulltext)
                % text(obj.schemPositions(i,1)+Xoff,obj.schemPositions(i,2)+Yoff, compNames{i})
                % if iscell(comps(i).paramExpressions) && numel(comps(i).paramExpressions) > 1
                %     text(obj.schemPositions(i,1)+Xoff,obj.schemPositions(i,2)+Yoff-10, strjoin(comps(i).paramExpressions( ~cellfun(@isempty,comps(i).paramExpressions)), '; ' ), ...
                %         'Color',[.5 .5 .5],'FontSize',10)
                % else
                %     text(obj.schemPositions(i,1)+Xoff,obj.schemPositions(i,2)+Yoff-10, comps(i).paramExpressions, ...
                %         'Color',[.5 .5 .5],'FontSize',10)
                % end
            end 
        end
    else
        warning('No data for schematic layout of components')

        f2 = figure('visible','off');

        %% block pulled from findNormalTreee -- should refactor it out
        nodes = unique([obj.origComponents(:).Nodes]);
        if ~isempty(find(strcmp(nodes,'0'),1))
            nodes = circshift(nodes, -1*(find(strcmp(nodes,'0'))-1)); % make '0' the first node
        end
        nodeMap = dictionary(nodes, 1:length(nodes));
        numNodes = reshape(nodeMap([obj.origComponents.Nodes]), [2,length(obj.origComponents),])';
    
        % Numerical nodes corresponding to components (not necessary here)
        typeMap = dictionary({'V','BV', 'MV','C','R','L','MI','BI','I', ...
            'E', 'F', 'Vm', 'Im', 'M', 'D'}, [1:9, 2, 8, 3, 7, 10:11]);
        typeOrder = typeMap({obj.origComponents.Type})';  
    
        subgraph = [typeOrder, numNodes , transpose(1:length(typeOrder))];
       
        origRowCount = size(subgraph,1);
        wireNum = max(typeOrder)+1;
        for i = origRowCount:-1:1
            if all(subgraph(i,:) == 0)
                subgraph(i,:) = [];
                continue;
            end
            newNodeNum = max(subgraph(:,2:3),[],'all')+1;
            origRow = subgraph(i,:);
            subgraph(i,:) = [];

            subgraph = [subgraph;
                wireNum origRow(2) newNodeNum Inf;
                origRow(1) newNodeNum newNodeNum+1 origRow(4);
                wireNum origRow(3) newNodeNum+1 Inf;];

            %Parallel elements
            parElements = all(subgraph(1:i-1,2:3) == origRow(2:3),2) | ...
                all(subgraph(1:i-1,2:3) == origRow(3:-1:2),2);
            if any(parElements)
                for j = find(parElements)
                    parRow = subgraph(j,:);
                    subgraph(j,:) = zeros(1,4);
                    subgraph = [subgraph;
                        wireNum newNodeNum newNodeNum+2 0;
                        parRow(1) newNodeNum+2 newNodeNum+3 parRow(4);
                        wireNum newNodeNum+3 newNodeNum+1 0];

                end
            end

                
        end



        G = graph(subgraph(:,2)', subgraph(:,3)', subgraph(:,4)');
        % G.Edges.Name = {obj.origComponents(G.Edges.Weight).Name}';
        Names = {};
        Names(~isinf(G.Edges.Weight) & G.Edges.Weight>0) = ...
            {obj.origComponents(G.Edges.Weight(~isinf(G.Edges.Weight) & G.Edges.Weight>0)).Name};
        Names(isinf(G.Edges.Weight)) = {'W'};
        Names(G.Edges.Weight ==0) = {'W'};
        compLocators = G.Edges.Weight;
        G.Edges.Weight(~isinf(G.Edges.Weight)  & G.Edges.Weight>0) = 1;
        G.Edges.Weight(~isinf(G.Edges.Weight)  & G.Edges.Weight==0) = 2;
        G.Edges.Weight(isinf(G.Edges.Weight)) = 3;
        G.Edges.Name = Names';


        
        %% 
    
        P = plot(G);
        layout(P,'force','WeightEffect','direct','UseGravity',false)
        [th,~] = cart2pol(P.XData(1), P.YData(1));
        rot = 3*pi/2-th;

        [ths,rs] = cart2pol(P.XData, P.YData);
        [P.XData, P.YData] = pol2cart(ths + rot, rs);
        P.YData(1) = min(P.YData)*1.25;

        nodeDists = (P.XData - P.XData').^2 + (P.YData - P.YData').^2;
        nodeDists = nodeDists + max(nodeDists,[],'all')*eye(numel(P.XData));
        minDist = min(nodeDists,[],'all');

        P.XData = P.XData/minDist;
        P.YData = P.YData/minDist;



        P.EdgeLabel = G.Edges.Name;

        

        for i = 1:numel(obj.origComponents)
            nodes = G.Edges.EndNodes(compLocators == i,:);
            compLoc = [mean(P.XData(nodes)) mean(P.YData(nodes))];
            compRot = round(mod(atan2(diff(P.YData(nodes)), diff(P.XData(nodes)))*180/pi,360)/90)*90;
            compSF = 100;

            obj.schemPositions(i,:) = [compLoc, compRot, compSF];

        end


        close(f2);

        plotSchematic(obj,fn)

        x=1
    end

    for i = 1:size(obj.schemWires,1)
        plot(obj.schemWires(i,[1 3]), -obj.schemWires(i,[2 4]),'k','LineWidth',1);
        plot(obj.schemWires(i,[1 3]), -obj.schemWires(i,[2 4]),'.k','LineWidth',1,'MarkerSize',10);
    end

    ax = f.Children(1);

    ax.XTick = [];
    ax.YTick = [];
    set(gca,'Visible','off')
end