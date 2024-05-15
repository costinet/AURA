function [A,B,C,D,I] = solveStateSpaceRepresentation(obj)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
%   Method adapted from L. O. Chua and P-M Lin, Compter Aided Analysis of
%   Electronic Circuits: Algorithms & Computational Techniques", Chapter 6 & 8.


    % Needed anymore?
    sources = obj.components(strcmp({obj.components.Type},'V') | strcmp({obj.components.Type},'I'));
    [~,IA,~] = intersect({obj.components.Name}, {sources.Name});
    obj.Element_Properties = [{obj.components.Name}' {obj.components.paramVals}'];
    obj.Element_Properties(IA,:) = [];
    obj.Component_Values = obj.Element_Properties;

    obj.findNormalTree();


    %% Below equivalent to old call to ABCD_num
    [almost_H] = obj.nodeloop_num(obj.NL,obj.NLnets);
    [H,s]=obj.hybridparse(almost_H,obj.SortedTree_cutloop,obj.SortedCoTree_cutloop);

    % Functions to find outputs
    [A,B,C,D,I,~,~,~,OutputNames,~,ConstantNames,OrderedNamesnum]=obj.loopfixAB_num(H,s,obj.NLnets,obj.SortedTree_cutloop,obj.SortedCoTree_cutloop);
    [C,D,~,~,StateNamesCD]=obj.loopfixCD_num(A,B,C,D,H,s,obj.NLnets,obj.SortedTree_cutloop,obj.SortedCoTree_cutloop);

    obj.StateNames = obj.NLnets(OrderedNamesnum,1);
    obj.OutputNamesCD = StateNamesCD;
    obj.OutputNames = OutputNames;
    obj.ConstantNames = ConstantNames;

    FETs = strcmp({obj.origComponents.Type} , 'M');
    diodes = strcmp({obj.origComponents.Type} , 'D');
    
    switches = FETs | diodes;

    obj.Switch_Names = {obj.origComponents(switches).Name};

end