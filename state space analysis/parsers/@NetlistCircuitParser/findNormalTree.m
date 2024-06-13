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
    nodeMap = dictionary(nodes, 1:length(nodes));
    numNodes = reshape(nodeMap([obj.components.Nodes]), [2,length(obj.components),])';

    %% Numerical nodes corresponding to components
    typeMap = dictionary({'V','BV', 'MV','C','R','L','MI','BI','I', ...
        'E', 'F', 'Vm', 'Im'}, [1:9, 2, 8, 3, 7]);
    typeOrder = typeMap({obj.components.Type})';  

    obj.NL = [typeOrder, numNodes , transpose(1:length(typeOrder))];
    obj.NLnets = [{obj.components.Name}' {obj.components.Nodes}' ];
    
    SortedRows = sortrows(obj.NL,1); % sorts rows in prefered tree order
    
    %% Find incidence matricies
    incidence = zeros(max(obj.NL(:,2:3),[],'all'),size(obj.NL,1)); % initialize incidence matrix
    
    idp = sub2ind(size(incidence), SortedRows(:,2), (1:length(SortedRows(:,2)))');
    idm = sub2ind(size(incidence), SortedRows(:,3), (1:length(SortedRows(:,3)))');
    
    incidence(idp)=1;
    incidence(idm)=-1;
    
    
    % Find the Reduced Incidence Matrix
    incidence_1 = incidence(2:end,:);
    
    % Reduce Incidence Matrix to Echelon Form
    REF = rref(incidence_1);

    if any(all(REF ==0,2))
        % Some rows are completely zero in the incidence matrix
        % Eliminating some nodes?  What does this mean?
        error(['Invalid circuit.  This may be becuase some elements are unconnected,' ...
            ' invalid components used, or some other error']);
        % incidence_1 = incidence_1(~all(REF ==0,2),:);
        % REF = rref(incidence_1);
    end
    
    [numRowsI,numColsI]=size(incidence_1);
    % Tree = 0;
    
    %% Find the Normal Tree and CoTree
    % Find the first one (1) in each row
    % The columns where these ones (1) occur are the braches within the normal
    % tree
    try
        colIndex = (REF==1).*repmat(1:numColsI,[numRowsI,1]);
        colIndex(colIndex <= 0) = inf;
        TreeRows = min(colIndex,[],2)';
    catch
        error('Possible Singular Incidence Matrix. It is possible that the circuit is not connected.')
    end
    
    
    % Find the remaining elements that are not in the normal tree
    % The branches are in the Cotree
    CoTreeRows = setdiff(1:numColsI,TreeRows);
    
    %% Orders Tree and CoTree indices correctly: E R G J
    
    % Adjust Branch Identification numbers from:
    % V BV MV  C  R  L MI BI  I
    % 1  2  3  4  5  6  7  8  9
    
    % to
    
    % Branch Identification numbers:
    % __________Tree__________  _________CoTree__________
    % E  E-B  E-M  E-C  E-L  R  G   J-C  J-L J-M  J-B  J
    % 1   2    3    4    5   6  7    8    9   10   11  12
    
    treeIDnumberMap = dictionary(1:9, [1:4 6 5 7:9]);
    coTreeIDnumberMap = dictionary(1:9, [1:3 8 7 9 10 11 12]);
   
    
    tree = [SortedRows(TreeRows(:),:), TreeRows(:)];
    tree(:,1) = treeIDnumberMap(tree(:,1));
    tree = sortrows(tree,1);
    
    coTree = [SortedRows(CoTreeRows(:),:), CoTreeRows(:)];
    coTree(:,1) = coTreeIDnumberMap(coTree(:,1));
    coTree = sortrows(coTree,1);
    
    assert(tree(end,1) <= 6, 'tree contains current sources, indicating the presence of a current net in the circuit.')
    assert(coTree(1,1) >= 7, 'CoTree contains voltage sources, indicating the presence of a voltage loop in the circuit.')
    
    %% Store resutls
    obj.Cutset = [REF(:,coTree(:,5))];
    obj.SortedTree_cutloop = tree;
    obj.SortedCoTree_cutloop = coTree;


end

