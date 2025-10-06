function [A,B,C,D,I] = solveStateSpaceRepresentation(obj)
%solveStateSpaceRepresentation solves state space matrices for
%circuitParser class.  
%   Requires that obj.components is populated.
%
%   Method adapted from L. O. Chua and P-M Lin, Compter Aided Analysis of
%   Electronic Circuits: Algorithms & Computational Techniques", Chapters 6 & 8.

    sources = obj.components(strcmp({obj.components.Type},'V') | strcmp({obj.components.Type},'I'));
    [~,IA,~] = intersect({obj.components.Name}, {sources.Name});
    obj.Element_Properties = [{obj.components.Name}' {obj.components.paramVals}'];
    obj.Element_Properties(IA,:) = [];
    obj.Component_Values = obj.Element_Properties;

    %% New functions
    obj.findNormalTree();
    [H,comps] = obj.hybridMatrix();
    [A,B,C,D,I,names] = obj.stateSpaceFromHybrid(H, comps);

    % new = [[H{1,1} H{3,1}; H{3,1} H{3,3}],[H{1,4}; H{3,4}]];


    %% New
    new = 1;
    if(new)
        if isempty(obj.StateNames)
            obj.StateNames = names{1};
        else
            assert(all(strcmp(obj.StateNames,names{1})),'When parsing a new switching subinterval, it appears that the states of the circuit are different than those from a previous subinterval')
        end
        obj.OutputNamesCD = names{2};
        obj.ConstantNames = names{3};
    else
    %% Old
        [almost_H] = obj.nodeloop_num();
        [Hold,s]=obj.hybridparse(almost_H,obj.SortedTree_cutloop,obj.SortedCoTree_cutloop);
        [Aold,Bold,Cold,Dold,Iold,~,~,~,OutputNames,~,ConstantNames,OrderedNamesnum]=obj.loopfixAB_num(Hold,s,obj.SortedTree_cutloop,obj.SortedCoTree_cutloop);
        [Cold,Dold,~,~,StateNamesCD]=obj.loopfixCD_num(Aold,Bold,Cold,Dold,Hold,s,obj.SortedTree_cutloop,obj.SortedCoTree_cutloop);

        A = Aold; B = Bold; C = Cold; D = Dold; I = Iold;
    
        obj.StateNames = {obj.components(OrderedNamesnum).Name}';
        obj.OutputNamesCD = StateNamesCD;
        obj.OutputNames = OutputNames;
        obj.ConstantNames = ConstantNames;
    end



    FETs = strcmp({obj.origComponents.Type} , 'M');
    diodes = strcmp({obj.origComponents.Type} , 'D');
    
    switches = FETs | diodes;

    obj.Switch_Names = {obj.origComponents(switches).Name};


    % H = obj.hybrid();
    % if ~all(H==almost_H,"all")
    %     error('New method broke')
    % end

end