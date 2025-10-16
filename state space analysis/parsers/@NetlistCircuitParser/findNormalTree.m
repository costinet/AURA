function findNormalTree(obj)
    %findNormalTree finds a normal tree spanning the circuit with a
    %priority towards voltage-source-like components in the tree
    %
    %   Method adapted from L. O. Chua and P-M Lin, Compter Aided Analysis of
    %   Electronic Circuits: Algorithms & Computational Techniques", Chapter 3.
    

    %% Interpret
    % Give nodes numeric identifiers
    nodes = unique([obj.components(:).Nodes]);
    if ~isempty(find(strcmp(nodes,'0'),1))
        nodes = circshift(nodes, -1*(find(strcmp(nodes,'0'))-1)); % make '0' the first node
    end
    nodeMap = dictionary(nodes, 1:length(nodes));
    numNodes = reshape(nodeMap([obj.components.Nodes]), [2,length(obj.components),])';

    %% Numerical nodes corresponding to components
    % typeMap = dictionary({'V','BV', 'MV','C','R','L','MI','BI','I', ...
    %     'E', 'F', 'Vm', 'Im'}, [1:9, 2, 8, 3, 7]);
    typeOrder = obj.typeMap({obj.components.Type})';  

    NL = [typeOrder, numNodes , transpose(1:length(typeOrder))];
    % obj.NLnets = [{obj.components.Name}' {obj.components.Nodes}' ];
 
    sortedComponents = sortrows(NL,1); % sorts components in prefered tree order
    
    %% Find incidence matrix
    incidence = zeros(numel(nodes),numel(obj.components));
    idp = sub2ind(size(incidence), sortedComponents(:,2), (1:length(sortedComponents(:,2)))');
    idm = sub2ind(size(incidence), sortedComponents(:,3), (1:length(sortedComponents(:,3)))');
    incidence(idp)=1;
    incidence(idm)=-1;   
    
    % Find the Reduced Incidence Matrix
    incidence = incidence(2:end,:);
    A = rref(incidence);

    if any(all(A ==0,2))
        % Some rows are completely zero in the incidence matrix
        if ~isempty(obj.netListDirectives)  && sum(all(A ==0,2)) == sum(strcmp({obj.netListDirectives.Type},'K'))
            % In transformer-isolated circuits, there will be unconnected
            % segments.  So this can continue, though the following may not
            % be the best wy to handle it.
            incidence = incidence(~all(A ==0,2),:);
            A = rref(incidence);
        else
            error(['Invalid circuit.  This may be becuase some elements are unconnected,' ...
                 ' invalid components used, or some other error']);
        end
    end
    
    [numRowsI,numColsI]=size(incidence);
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
    % % E  E-B  E-M  E-C  R  E-L  J-C  G   J-L J-M  J-B  J
    % % 1   2    3    4   5   6    7   8    9   10   11  12

    
    tree = [sortedComponents(treeRows(:),:), treeRows(:)];
    tree(:,1) = obj.treeIDnumberMap(tree(:,1));
    tree = sortrows(tree,1);
    
    coTree = [sortedComponents(coTreeRows(:),:), coTreeRows(:)];
    coTree(:,1) = obj.coTreeIDnumberMap(coTree(:,1));
    coTree = sortrows(coTree,1);
    
    assert(tree(end,1) <= 6, 'tree contains current sources, indicating the presence of a current net in the circuit.')
    assert(coTree(1,1) >= 7, 'coTree contains voltage sources, indicating the presence of a voltage loop in the circuit.')

    if(0)
        % These are inaccessible.  These loops/cutsets are now handled,
        % including state-source dependence, so the warning may not be
        % necessary (though could still be informative).
        if any(tree(:,1) == 5)
            warning(['L-I net present.  ' strjoin({obj.components(coTree(tree(:,1) == 5,4)).Name}, ', ') ' will be dependent state(s)'])
        end
        if any(coTree(:,1) == 8)
            warning(['C-V loop present.  ' strjoin({obj.components(coTree(coTree(:,1) == 8,4)).Name}, ', ') ' will be dependent state(s)']);
        end
    end
    
    %% Store results
    AL = A(:,coTree(:,5));

    AT = A(:,sort(tree(:,5)));
    assert(all(AT == eye(size(AT)),'all'), 'Tree incorrect')

    obj.AL = AL;
    obj.tree = tree;
    obj.coTree = coTree;

    return

    %% Unreachable, just for debugging:
    
    A = [AT, AL];
    BT = (-AT\AL)';
    B = [BT, eye(size(BT,1))];
    D = AT\A;

    rowMap = [treeRows, coTreeRows];

    % fundamental loops:
    for i = 1:size(B,1)
        compLoc = rowMap(B(i,:)~=0);
        {obj.components(sortedComponents(compLoc,end)).Name}
        loopNodes = [[obj.components(sortedComponents(compLoc,end)).Nodes]'];
        assert(length(loopNodes) == 2*length(unique(loopNodes)), 'This is not a loop');
    end

    % fundamental cutsets:
    for j = 1:size(D,1)
        compLoc = rowMap(D(j,:)~=0);
        {obj.components(sortedComponents(compLoc,end)).Name}
    end
end

