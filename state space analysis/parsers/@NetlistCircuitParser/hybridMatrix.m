function [H,remainingComps] = hybridMatrix(obj)
%hybridMatrix solves hybrid matrix for resistve N-port, including
%dependent sources
%   Detailed explanation goes here
    

    DL = obj.AL;
    tree = obj.tree;
    coTree = obj.coTree;

    treeResistors = tree(:,1) == obj.treeIDnumberMap(obj.typeMap({'R'}));
    coTreeResistors = coTree(:,1) == obj.coTreeIDnumberMap(obj.typeMap({'R'}));

    Rvals = [obj.components(tree(treeResistors,4)).paramVals];
    Gvals = [obj.components(coTree(coTreeResistors,4)).paramVals];
    ZR = diag(Rvals);
    ZG = diag(Gvals);
    YR = diag(1./Rvals);
    YG = diag(1./Gvals);

    EG = {~treeResistors, coTreeResistors};
    EJ = {~treeResistors, ~coTreeResistors};
    RG = {treeResistors, coTreeResistors};
    RJ = {treeResistors, ~coTreeResistors};

    %% Solve Z, Y, and H
    % See Eqs (6-19) & (6-20) in Chua [ISBN 0131654152]
    % This is for the purely-resistive m-port

    Z = ZG+(DL(RG{:}).')*ZR*DL(RG{:});
    Y = (YR)+DL(RG{:})*(YG)*(DL(RG{:}).');
    H = [-DL(EG{:})*(Z\(DL(EG{:})')),... %H_EE
         DL(EG{:})*(Z\DL(RG{:})')*ZR*DL(RJ{:})-DL(EJ{:}); ... %H_EJ
         (DL(EJ{:}).')-((DL(RJ{:})')*(Y\DL(RG{:}))*(YG)*(DL(EG{:})')), ... %H_JE
         -(DL(RJ{:})')*(Y\DL(RJ{:}))  ]; %H_JJ

    % Hcomponents = {obj.components(tree(~treeResistors,4)).Name, '||', obj.components(coTree(~coTreeResistors,4)).Name};

    remainingComps = [obj.components(tree(~treeResistors,4)), obj.components(coTree(~coTreeResistors,4))];
    remCompsTypes = {remainingComps.Type};

    %%
    states = strcmp(remCompsTypes,'C') | strcmp(remCompsTypes,'L');
    measSources = strcmp(remCompsTypes,'Vm') | strcmp(remCompsTypes,'Im');
    depSources = strcmp(remCompsTypes,'E') | strcmp(remCompsTypes,'F') | strcmp(remCompsTypes,'G') | strcmp(remCompsTypes,'H');
    indepSources = strcmp(remCompsTypes,'V') | strcmp(remCompsTypes,'I');

    assert(all(states + measSources + depSources + indepSources == 1), 'Unknown component in netlist');

    % states = states | measSources; % | depSources; %Everything we want to solve for, except indep sources (will be in s).

    H = {H(states,states), H(states,depSources), H(states,measSources), H(states,indepSources);
        H(depSources,states), H(depSources,depSources), H(depSources,measSources), H(depSources,indepSources);
        H(measSources,states), H(measSources,depSources), H(measSources,measSources), H(measSources,indepSources);
        H(indepSources,states), H(indepSources,depSources), H(indepSources,measSources), H(indepSources,indepSources) };
    
    
        
    %% Find coefficients of dependent sources 
    K = zeros(sum(depSources), sum(measSources));
    depSourceComps = remainingComps(depSources);
    for i=1:numel(depSourceComps)
        comp = depSourceComps(i);
        if numel({comp.paramVals{1:2:end}}) == 1
            refSource = strcmp({remainingComps.Name}, comp.paramVals{1});
            assert(any(measSources & refSource), 'Currently, A dependent source must only be dependent on the voltage/current of a measurement source')
            refIndex = find(find(measSources)==find(refSource));
            K(i,refIndex) = comp.paramVals{2};
        else
            error('This has not been updated, yet, after being copied from hybridparse.')
            %We don't currently allow any components corresponding to this
            %else condition, but the old code is left in case this comes up
            %later.  It has not been updated to work with the current
            %class.
            % for j = 1:numel({obj.components(depSources(i)).paramVals{1:2:end}})
            %     refSource = find(strcmp({obj.components.Name}, obj.components(depSources(i)).paramVals{1 + 2*(j-1)}));
            %     assert(~isempty(intersect(depSources,refSource)), 'Currently, only linked dependent sources are possible.  A dependent source must only be dependent on the voltage/current of another dependent source')
            %     [~,refIndex] = intersect(depSources,refSource);
            %     K(i,refIndex) = obj.components(depSources(i)).paramVals{2*j};
            % end
        end
    end
    
    if ~isempty(K)
        F = K*H{3,2};

        if cond(eye(size(F))-F) > 1e8
            warning('Singular Matrix when trying to solve for dependent sources. This is likely due to poor port assignment.  Try changing the ordering of winding inductors in the K (coupling) directive.')
        end
    end

    %% Eliminate Dependent Sources
    % See Eqs (6-49) & (6-50) in Chua [ISBN 0131654152]

    I22 = eye(size(H{3,2}*K));
    H{1,1} = H{1,1} + H{1,2}*K*((I22-H{3,2}*K)\H{3,1}); %(H)
    H{1,4} = H{1,4} + H{1,2}*K*((I22-H{3,2}*K)\H{3,4}); %(s)
    
    H{3,1} = H{3,1} + H{3,2}*K*((I22-H{3,2}*K)\H{3,1}); 
    H{3,4} = H{3,4} + H{3,2}*K*((I22-H{3,2}*K)\H{3,4}); 

    if any(isnan(H{1,1}),'all')
        error('Singular Matrix when trying to solve for dependent sources. This is likely due to poor port assignment.  Try changing the ordering of winding inductors in the K (coupling) directive.')
    end

    remainingComps = {remainingComps(states), remainingComps(depSources), remainingComps(measSources), remainingComps(indepSources)};


end