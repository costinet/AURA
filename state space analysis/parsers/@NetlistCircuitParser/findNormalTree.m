function findNormalTree(obj)
    %findNormalTree 
    %   Detailed explanation goes here
    
    %note: replaces cutset_loop_num

    %% Interpret
    % Give nodes numeric identifiers
    nodes = unique([obj.components(:).Nodes]);
    if ~isempty(find(strcmp(nodes,'0'),1))
        nodes = circshift(nodes, -1*(find(strcmp(nodes,'0'))-1)); % make '0' the first node
    end
    obj.nodeMap = dictionary(nodes, 1:length(nodes));

    % [nodeMap, shorts] = findShorts(obj, nodeMap);
    numNodes = reshape(obj.nodeMap([obj.components.Nodes]), [2,length(obj.components),])';

    %% Numerical nodes corresponding to components
    % typeMap = dictionary({'V','BV', 'MV','C','R','L','MI','BI','I', ...
    %     'E', 'F', 'Vm', 'Im'}, [1:9, 2, 8, 3, 7]);
    typeOrder = obj.typeMap({obj.components.Type})';  

    obj.NL = [typeOrder, numNodes , transpose(1:length(typeOrder))];
    % obj.NLnets = [{obj.components.Name}' {obj.components.Nodes}' ];
 
    sortedComponents = sortrows(obj.NL,1); % sorts components in prefered tree order
    
    %% Find incidence matricies
    % nodeCount = numel(unique(obj.NL(:,2:3)));
    % incidence = zeros(nodeCount,size(obj.NL,1)); % initialize incidence matrix
    % incidence = zeros(max(obj.NL(:,2:3),[],'all'),size(obj.NL,1)); % initialize incidence matrix
    incidence = zeros(numel(nodes),numel(obj.components));
    
    
    idp = sub2ind(size(incidence), sortedComponents(:,2), (1:length(sortedComponents(:,2)))');
    idm = sub2ind(size(incidence), sortedComponents(:,3), (1:length(sortedComponents(:,3)))');
    
    incidence(idp)=1;
    incidence(idm)=-1;

    % Remove rows corresponding to shorted elements
    % incidence = incidence(sort(unique(nodeMap.values)),:);
    
    
    % Find the Reduced Incidence Matrix
    incidence_1 = incidence(2:end,:);
    
    % Reduce Incidence Matrix to Echelon Form
    A = rref(incidence_1);

    if any(all(A ==0,2))
        % Some rows are completely zero in the incidence matrix
        if ~isempty(obj.netListDirectives)  && sum(all(A ==0,2)) == sum(strcmp({obj.netListDirectives.Type},'K'))
            % In transformer-isolated circuits, there will be unconnected
            % segments.  So this can continue, though the following may not
            % be the best wy to handle it.
            incidence_1 = incidence_1(~all(A ==0,2),:);
            A = rref(incidence_1);
        else
            error(['Invalid circuit.  This may be becuase some elements are unconnected,' ...
                 ' invalid components used, or some other error']);
        end
    end
    
    [numRowsI,numColsI]=size(incidence_1);
    % Tree = 0;
    
    %% Find the Normal Tree and CoTree
    % Find the first one (1) in each row
    % The columns where these ones (1) occur are the braches within the normal
    % tree
    try
        colIndex = (A==1).*repmat(1:numColsI,[numRowsI,1]);
        colIndex(colIndex <= 0) = inf;
        treeRows = min(colIndex,[],2)';
    catch
        error('Possible Singular Incidence Matrix. It is possible that the circuit is not connected.')
    end
    
    % Find the remaining elements that are not in the normal tree
    % The branches are in the Cotree
    coTreeRows = setdiff(1:numColsI,treeRows);


    
    % %% Orders Tree and CoTree indices correctly: E R G J
    % 
    % % Adjust Branch Identification numbers from:
    % % V BV MV  C  R  L MI BI  I
    % % 1  2  3  4  5  6  7  8  9
    % 
    % % to
    % 
    % % Branch Identification numbers:
    % % __________Tree__________  _________CoTree__________
    % % E  E-B  E-M  E-C  E-L  R  G   J-C  J-L J-M  J-B  J
    % % 1   2    3    4    5   6  7    8    9   10   11  12
    % 
    % treeIDnumberMap = dictionary(1:9, [1:4 6 5 7:9]);
    % coTreeIDnumberMap = dictionary(1:9, [1:3 8 7 9 10 11 12]);
   
    
    tree = [sortedComponents(treeRows(:),:), treeRows(:)];
    tree(:,1) = obj.treeIDnumberMap(tree(:,1));
    tree = sortrows(tree,1);
    
    coTree = [sortedComponents(coTreeRows(:),:), coTreeRows(:)];
    coTree(:,1) = obj.coTreeIDnumberMap(coTree(:,1));
    coTree = sortrows(coTree,1);
    
    assert(tree(end,1) <= 6, 'tree contains current sources, indicating the presence of a current net in the circuit.')
    assert(coTree(1,1) >= 7, 'CoTree contains voltage sources, indicating the presence of a voltage loop in the circuit.')

    if(0)
        % These warnings are suppressed so they don't show up repeatedly on
        % each switching interval parse
        if any(tree(:,1) == 5)
            warning(['L-I net present.  ' strjoin({obj.components(coTree(tree(:,1) == 5,4)).Name}, ', ') ' will be dependent state(s)'])
        end
        if any(coTree(:,1) == 8)
            warning(['C-V loop present.  ' strjoin({obj.components(coTree(coTree(:,1) == 8,4)).Name}, ', ') ' will be dependent state(s)']);
        end
    end
    
    %% Store results
    AT = A(:,tree(:,5));
    AL = A(:,coTree(:,5));
    obj.Cutset = AL;
    obj.SortedTree_cutloop = tree;
    obj.SortedCoTree_cutloop = coTree;

    return

    % Not evaluated, just for debugging:
    A = [AT, AL];
    BT = (-AT\AL)';
    B = [BT, eye(size(BT,1))];
    D = AT\A;

    % fundamental loops:
    for i = 1:size(A,1)
        obj.components(sortedComponents(B(i,:)~=0,end)).Name
    end

    % fundamental cutsets:
    for j = 1:size(D,1)
        {obj.components(sortedComponents(D(j,:)~=0,end)).Name}
    end
end

