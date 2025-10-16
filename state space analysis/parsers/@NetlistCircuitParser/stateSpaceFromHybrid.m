function [A,B,C,D,I,names] = stateSpaceFromHybrid(obj, H, comps)
%stateSpaceFromHybrid solves state space matrices from hybrid matrics
%   [A,B,C,D,I,names] = stateSpaceFromHybrid(obj, H, comps)
%
%   Method adapted from L. O. Chua and P-M Lin, Compter Aided Analysis of
%   Electronic Circuits: Algorithms & Computational Techniques", Chapter 8.

    tree = obj.tree;
    coTree =  obj.coTree;

    M = diag([comps{1}.paramVals]);

    dependentStates = obj.components([tree(tree(:,1) == obj.treeIDnumberMap(obj.typeMap({'L'})),4);  ...
        coTree(coTree(:,1) == obj.coTreeIDnumberMap(obj.typeMap({'C'})),4)]);

    origIndex = zeros(numel(comps{1}),1);
    for i = 1:numel(comps{1})
        origIndex(i) = find(strcmp(comps{1}(i).Name, {obj.components.Name}));
    end
    if ~isempty(dependentStates)

        fvec = false(1,length(H{1,1}));
        [~,Itree] = intersect(origIndex, tree(:,4));
        [~,IcoTree] = intersect(origIndex, coTree(:,4));
        
        treecomps = fvec; treecomps(Itree) = true;
        coTreecomps = fvec; coTreecomps(IcoTree) = true;

        treeCaps = strcmp({comps{1}.Type},'C') & treecomps;
        coTreeCaps = strcmp({comps{1}.Type},'C') & coTreecomps;

        coTreeInductors = strcmp({comps{1}.Type},'L') & coTreecomps;
        treeInductors = strcmp({comps{1}.Type},'L') & treecomps;

        % F = {H{4,4}, H{4,1}(:,coTreeInductors), [], H{4,1}(:,coTreeCaps);
        %     H{1,4}(treeCaps,:), H{1,1}(treeCaps,coTreeInductors), [], H{1,1}(treeCaps,coTreeCaps);
        %     [], [], [], [];
        %     H{1,4}(treeInductors,:), H{1,1}(treeInductors,coTreeInductors), [], H{1,1}(treeInductors,coTreeCaps)      
        % };

        
        F = {H{1,1}(treeCaps,treeCaps), H{1,1}(treeCaps,coTreeCaps), H{1,1}(treeCaps,coTreeInductors), H{1,1}(treeCaps,treeInductors); ...
            H{1,1}(coTreeCaps,treeCaps), H{1,1}(coTreeCaps,coTreeCaps), H{1,1}(coTreeCaps,coTreeInductors), H{1,1}(coTreeCaps,treeInductors); ...
            H{1,1}(coTreeInductors,treeCaps), H{1,1}(coTreeInductors,coTreeCaps), H{1,1}(coTreeInductors,coTreeInductors), H{1,1}(coTreeInductors,treeInductors); ...
            H{1,1}(treeInductors,treeCaps), H{1,1}(treeInductors,coTreeCaps), H{1,1}(treeInductors,coTreeInductors), H{1,1}(treeInductors,treeInductors)       
        };
        

        CJ = diag([comps{1}(treeCaps).paramVals]);
        CL = diag([comps{1}(coTreeCaps).paramVals]);
        LJ = diag([comps{1}(treeInductors).paramVals]);
        LL = diag([comps{1}(coTreeInductors).paramVals]);

        rc = [sum(treeCaps), sum(coTreeCaps), sum(coTreeInductors), sum(treeInductors)];

        % M0 = [CJ + F{2,4}*CL*(F{2,4}'), zeros(sum(treeCaps),sum(coTreeInductors));
        %     zeros(sum(coTreeInductors), sum(treeCaps)), LL + (F{4,2}')*LJ*F{4,2}];
        % A0 = [-F{2,2}; F{2,2}']
        

        % M0 = [eye(rc(1)), -F{1,2}, zeros(rc(1),rc(3)), -F{1,4};
        %      zeros(rc(2), sum(rc)); % coTree Caps
        %      zeros(rc(3),rc(1)), -F{3,2}, eye(rc(3)), -F{3,4};
        %      zeros(rc(4), sum(rc))]; %tree Inductors
        % A0 = [F{1,1}, zeros(rc(1),rc(2)), F{1,3}, zeros(rc(1),rc(4));
        %     F{2,1}, -eye(rc(2)), zeros(rc(2),rc(3)), zeros(rc(2),rc(4));
        %     F{3,1}, zeros(rc(3),rc(2)), F{3,3}, zeros(rc(3),rc(4));
        %     zeros(rc(4),rc(1)), zeros(rc(4),rc(2)), F{4,3}, -eye(rc(4),rc(4))];
        % B0 = H{1,4};

        % rref([M0,A0,B0])

        % By the way these are forumlated, we can eliminate the rows from
        % the dependent components.  The dependent cap voltages wont' show
        % up anywhere and the dependent inductor currents won't be
        % anywhere.  But, the cap currents and inductor voltages are still
        % in M0.  The algebraic equation needs to be differentiated and
        % plugged in to eliminate them from M0

        % M0 = [eye(rc(1)), -F{1,2}, zeros(rc(1),rc(3)), -F{1,4};
        %      zeros(rc(3),rc(1)), -F{3,2}, eye(rc(3)), -F{3,4}];
        % A0 = [F{1,1}, zeros(rc(1),rc(2)), F{1,3}, zeros(rc(1),rc(4));
        %     F{3,1}, zeros(rc(3),rc(2)), F{3,3}, zeros(rc(3),rc(4))];
        % B0 = H{1,4}(treeCaps | coTreeInductors, :);

        

        % M0 = [eye(rc(1)), -F{1,2}, zeros(rc(1),rc(3)), -F{1,4};
        %      zeros(rc(3),rc(1)), -F{3,2}, eye(rc(3)), -F{3,4}];
        % A0 = [F{1,1}, zeros(rc(1),rc(2)), F{1,3}, zeros(rc(1),rc(4));
        %     F{3,1}, zeros(rc(3),rc(2)), F{3,3}, zeros(rc(3),rc(4))];
        % B0 = H{1,4}(treeCaps | coTreeInductors, :);

        % [F{2,1}, -eye(rc(2))] H{1,4}(coTreeCaps,:)

        % M0Try1 =  [(eye(rc(1)) -F{1,2}*F{2,1}),  F{1,4}*F{4,3};
        %        -F{3,2}*F{2,1}, (eye(rc(3)) - F{3,4}*F{4,3})];


        ICJ = diag(1./[comps{1}(treeCaps).paramVals]);
        ILL = diag(1./[comps{1}(coTreeInductors).paramVals]);

        M0 =  [(eye(rc(1)) -F{1,2}*CL*F{2,1}*ICJ),  -F{1,4}*LJ*F{4,3}*ILL;
               -F{3,2}*CL*F{2,1}*ICJ, (eye(rc(3)) - F{3,4}*LJ*F{4,3}*ILL)];
        A0 = [F{1,1}, F{1,3};
            F{3,1},F{3,3}];
        B0 = H{1,4}(treeCaps | coTreeInductors, :) + 0; % differentiation eliminates impact of indep sources.

        M = diag([comps{1}(treeCaps | coTreeInductors).paramVals]);

        % A = M\(M0\A0);
        % B = M\(M0\B0);

        %% Add back in dependent States
        % nDepSt = sum(coTreeCaps) + sum(treeInductors);

        % I = eye(size(A));
        % I = [I, zeros(size(I,1),nDepSt);
        %         F{2,1}, zeros(rc(2),rc(3)), zeros(rc(2)), zeros(rc(2),rc(4));
        %         zeros(rc(4),rc(1)), F{4,3}, zeros(rc(4),rc(2)), zeros(rc(4))];

        %This was incorrect because treeInductors actually come *first*
        % I = [eye(rc(1)), zeros(rc(1), sum(rc(2:end)));
        %         F{2,1}, zeros(rc(2),rc(3)), zeros(rc(2)), zeros(rc(2),rc(4));
        %         zeros(rc(3),sum(rc(1:2))), eye(rc(3)), zeros(rc(3),rc(4));
        %         zeros(rc(4),rc(1)), F{4,3}, zeros(rc(4),rc(2)), zeros(rc(4))];

        I = double(diag(treeCaps | coTreeInductors));
        I(coTreeCaps,treeCaps) = F{2,1};
        I(treeInductors,coTreeInductors) = F{4,3};

        % If there is a C-V loop involving a V, we need a I-like matrix for
        % the sources as well
        % BI = zeros(size(B));
        % BI = [BI;  H{1,4}(coTreeCaps | treeInductors, :)];
        BI = zeros(size(I,1),size(H{4,4},1));
        BI(coTreeCaps | treeInductors, :) =  H{1,4}(coTreeCaps | treeInductors, :);

        I = [I BI];

        if obj.useForcedDependentStates
            % This uses a fictitious resistor to force the
            % dependent states to their algebraic equations.  It should
            % generally lead to nicer plotted waveforms, but does add in
            % non-real dynamics to the system model

            rb = .001;
            % A = [A, zeros(size(A,1),nDepSt);
            %     F{2,1}/rb, zeros(rc(2),rc(3)), -eye(rc(2))/rb, zeros(rc(2),rc(4));
            %     zeros(rc(4),rc(1)), F{4,3}, zeros(rc(4),rc(2)), -eye(rc(4))/rb];
            % B = [B;  H{1,4}(coTreeCaps | treeInductors, :)/rb];

             % A = [A(1:rc(1),1:rc(1)), zeros(rc(1),rc(2)), A(1:rc(1),rc(1)+1:rc(1)+rc(3)), zeros(rc(1), rc(4));
             %    F{2,1}/rb, zeros(rc(2),rc(3)), -eye(rc(2))/rb, zeros(rc(2),rc(4));
             %    A(rc(1)+1:rc(1)+rc(3),1:rc(1)), zeros(rc(3),rc(2)), A(rc(1)+1:rc(1)+rc(3),rc(1)+1:rc(1)+rc(3)), zeros(rc(3), rc(4));
             %    zeros(rc(4),rc(1)), F{4,3}, zeros(rc(4),rc(2)), -eye(rc(4))/rb];
             A = zeros(sum(rc));
             A(treeCaps | coTreeInductors, treeCaps | coTreeInductors) = M\(M0\A0);
             A(coTreeCaps | treeInductors,:) = [
                 F{2,1}/rb, zeros(rc(2),rc(3)), -eye(rc(2))/rb, zeros(rc(2),rc(4))
                 zeros(rc(4),rc(1)), F{4,3}, zeros(rc(4),rc(2)), -eye(rc(4))/rb
                 ];
             B = zeros(sum(rc),size(H{4,4},1));
             B(treeCaps | coTreeInductors,:) = M\(M0\B0);
             B(coTreeCaps | treeInductors,:) = diag([diag(CL)', diag(LJ)'])\H{1,4}(coTreeCaps | treeInductors, :)/rb;

        else
            % This just copies the derivatives of the dependent states to match the
            % derivatives of other states according to the algebraic
            % constraint equations. This is correct dynamics, but the DC
            % operating point is unconstrained, so they may show up on plots
            % with random bias. 

            % A = [A, zeros(size(A,1),nDepSt);
            %     F{2,1}*A(find(treeCaps), find(treeCaps)), zeros(rc(2),rc(3)), zeros(rc(2)), zeros(rc(2),rc(4));
            %     zeros(rc(4),rc(1)), F{4,3}*A(find(coTreeInductors)-rc(2), find(coTreeInductors)-rc(1)), zeros(rc(4),rc(2)), zeros(rc(4))];
            % B = [B;  diag([diag(CL)', diag(LJ)'])\H{1,4}(coTreeCaps | treeInductors, :)];

             A = zeros(sum(rc));
             A(treeCaps | coTreeInductors, treeCaps | coTreeInductors) = M\(M0\A0);
             A(coTreeCaps | treeInductors,:) = [
                 % F{2,1}*A(find(treeCaps), find(treeCaps)), zeros(rc(2),rc(3)), -eye(rc(2)), zeros(rc(2),rc(4)) %I thought the find was necessary due to a size mismatch, but may be wrong
                 F{2,1}*A(treeCaps, treeCaps), zeros(rc(2),rc(3)), -eye(rc(2)), zeros(rc(2),rc(4)) 
                 zeros(rc(4),rc(1)), F{4,3}*A(find(coTreeInductors)-rc(2), find(coTreeInductors)-rc(1)), zeros(rc(4),rc(2)), -eye(rc(4))
                 ];
             B = zeros(sum(rc),size(H{4,4},1));
             B(treeCaps | coTreeInductors,:) = M\(M0\B0);
             B(coTreeCaps | treeInductors,:) = diag([diag(CL)', diag(LJ)'])\H{1,4}(coTreeCaps | treeInductors, :);

        end

        %% Output Equations
        C0 = H{3,1};
        C = zeros(size(C0));
        C(:, treeCaps | coTreeInductors) = C0(:, treeCaps | coTreeInductors);
        
        % C = C + H{3,1}(:, coTreeCaps)*(CL*F{2,1})*A(treeCaps,:)...
        %     + H{3,1}(:, treeInductors)*LJ*F{4,3}*A(coTreeInductors,:);

        C = C + C0(:, coTreeCaps)*(CL*F{2,1})*A(treeCaps,:)...
            + C0(:, treeInductors)*LJ*F{4,3}*A(coTreeInductors,:);
        
        D0 = H{3,4};
        % D = D + H{3,4}(coTreeCaps,:)*(CL*F{2,1})*B(treeCaps,:) + ...
        %     H{3,4}(treeInductors,:)*LJ*F{4,3}*B(coTreeInductors,:); %% Problem Here.

        % D = zeros(size(D0));
        % D(treeCaps | coTreeInductors,:)  = D0(treeCaps | coTreeInductors,:);
        D = D0 + C0(:, coTreeCaps)*(CL*F{2,1})*B(treeCaps,:) + ...
            C0(:, treeInductors)*LJ*F{4,3}*B(coTreeInductors,:); 

        % names = {{comps{1}(treeCaps | coTreeInductors).Name, comps{1}(~treeCaps & ~coTreeInductors).Name},...
        %     {comps{3}.Name}, {comps{4}.Name}};
        names = {{comps{1}.Name}, {comps{3}.Name}, {comps{4}.Name}};
        

    else

        A = M\H{1,1};
        B = M\H{1,4};
        C = H{3,1};
        D = H{3,4};
        I = eye(size(H{1,1}));

        names = {{comps{1}.Name}, {comps{3}.Name}, {comps{4}.Name}};
    end
    
end