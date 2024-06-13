function [H, tree, coTree, nNL] = hybrid(obj)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    %% Numerical nodes corresponding to components
    nodes = unique([obj.components(:).Nodes]);
    if ~isempty(find(strcmp(nodes,'0'),1))
        nodes = circshift(nodes, -1*(find(strcmp(nodes,'0'))-1)); % make '0' the first node
    end
    nodeMap = dictionary(nodes, 1:length(nodes));
    numNodes = reshape(nodeMap([obj.components.Nodes]), [2,length(obj.components),])';

    %% Preferred order for elements in tree
    typeMap = dictionary({'V','BV', 'MV','C','R','L','MI','BI','I', ...
        'E', 'F', 'Vm', 'Im'}, [1:9, 2, 8, 3, 7]);
    mapType = dictionary(1:9,{'V','BV', 'MV','C','R','L','MI','BI','I'});
    typeOrder = typeMap({obj.components.Type})';  

    %% Numerical Netlist
    % [order, netConnections(1:2), originalPosition]
    nNL = [typeOrder, numNodes , transpose(1:numel(obj.components))];
    nNL = sortrows(nNL,1); % sorts rows in prefered tree order

    %% Incidence matrix
    I = zeros(length(nodes),numel(obj.components));
    I(sub2ind(size(I),nNL(:,2)',1:size(nNL,1))) = 1;
    I(sub2ind(size(I),nNL(:,3)',1:size(nNL,1))) = -1;

    % reduce, and eliminate one node (nominally ground)
    rI = rref(I(2:end,:));
    % rAdj = 

    %% Debugging visualizations:
    % disInc(I,obj.components(nNL(:,end)),nodes)
    % disInc(rI,obj.components(nNL(:,end)),nodes(2:end))

    assert(~any(all(rI ==0,2)), ['Invalid circuit.  This may be becuase some elements are unconnected,' ...
            ' invalid components used, or some other error']);

    %% Find the Normal Tree and CoTree
    % Find the first one (1) in each row
    % The columns where these ones (1) occur are the braches within the normal
    % tree
    try
        colIndex = (rI==1).*repmat(1:size(rI,2),[size(rI,1),1]);
        colIndex(colIndex <= 0) = inf;
        TreeRows = min(colIndex,[],2)';
    catch
        error('Possible Singular Incidence Matrix. It is possible that the circuit is not connected.')
    end

    % Find the remaining elements that are not in the normal tree
    % The branches are in the Cotree
    CoTreeRows = setdiff(1:size(rI,2),TreeRows);

    tree = nNL(TreeRows,:);
    coTree = nNL(CoTreeRows,:);
    cutSet = rI(:,CoTreeRows);

    assert(all(tree(:,1)<=6 ), 'tree contains current sources, indicating the presence of a current net in the circuit.');
    assert(all(coTree(:,1)>=4), 'CoTree contains voltage sources, indicating the presence of a voltage loop in the circuit.');
    % NOT Checked: are there any C-V loops?  L-I nets?

    %% 
    % Find all resistors in the tree
    Rlocs = tree(:,1)==typeMap({'R'});
    Z_R = diag([obj.components(tree(Rlocs,end)).paramVals]);
    Y_R =  diag(1./[obj.components(tree(Rlocs,end)).paramVals]);

    % Find all resistors in the coTree
    Glocs = coTree(:,1)==typeMap({'R'});
    Z_G = diag([obj.components(coTree(Glocs,end)).paramVals]);
    Y_G =  diag(1./[obj.components(coTree(Glocs,end)).paramVals]);

    %% Parsing D_L Matrix
    D_EG = cutSet(~Rlocs,Glocs);
    D_EJ = cutSet(~Rlocs,~Glocs);
    D_RG = cutSet(Rlocs,Glocs);
    D_RJ = cutSet(Rlocs,~Glocs);

    if isempty(D_EG) || isempty(D_EJ) || isempty(D_RG) || isempty(D_RJ)
        x=1;
    end

    %% Solve for Z and Y Matrices
    Z = Z_G+(D_RG.')*Z_R*D_RG;
    Y = Y_R+D_RG*Y_G*(D_RG.');

    H_EE = -(D_EG*(Z\(D_EG.')));
    H_EJ = D_EG*(Z\D_RG')*Z_R*D_RJ-D_EJ;
    H_JE = (D_EJ.')-((D_RJ')*(Y\D_RG)*Y_G*(D_EG'));
    H_JJ = -(D_RJ')*(Y\D_RJ);

    H=[H_EE,H_EJ;H_JE,H_JJ];
   return

end

    



function disInc(I,components,nodes)
    T = array2table(I,'RowNames',nodes,'VariableNames',{components.Name});
    disp(T)
end