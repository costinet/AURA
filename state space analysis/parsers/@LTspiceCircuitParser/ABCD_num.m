function [A,B,C,D,I] = ABCD_num(obj,Switch_Resistors,SW,Switch_Sequence)
% ABCD_num creates takes a NETlist file from LTSpice and creates the associated
% ABCD matrices or the associated matrices need to get ABCD if a symbolic
% formulation cannot be solved in a descent amount of time.
% [A,B,C,D] = ABCD_num(top.Switch_Resistors,top.Switch_Resistor_Values,top.Switch_Sequence)
% all inputs found in topology class
% A NetListParse class is required with the associated filename defined
%
% filename contains a string of the netlist filename from LTSpice Example:
% boost.net


%     _   _   _  ____    _    
%    / \ | | | |/ _  |  / \   
%   / _ \| | | | (_| | / _ \  
%  / ___ | |_| |> _  |/ ___ \ 
% /_/   \_\___//_/ |_/_/   \_\

% for i = 1:length(Switch_Sequence)
% Switch_Value{i} = SW(Switch_Sequence(i)+1,i);
% end

Switch_Value = diag(SW(1:length(Switch_Sequence),Switch_Sequence+1));
Switch_Value = num2cell(Switch_Value);

% obj.Component_Values(end-length(Switch_Resistors)+1:end,:) = [Switch_Resistors Switch_Value'];

%% Is this the whole reason why he needs them at the end?
[~,IA,IB] = intersect(obj.Component_Values(:,1), Switch_Resistors);
obj.Component_Values(IA,2) = Switch_Value(IB)';
obj.Component_Values = [obj.Component_Values; obj.Component_Values(IA,:)];
obj.Component_Values(IA,:) = [];



i = 1;
    % numeric solution
    [almost_H] = obj.nodeloop_num(obj.NL,obj.NLnets);
    [H,s]=obj.hybridparse(almost_H,obj.SortedTree_cutloop,obj.SortedCoTree_cutloop);

    % Functions to find outputs
    [A,B,C,D,I,~,~,~,OutputNames,~,ConstantNames,OrderedNamesnum]=obj.loopfixAB_num(H,s,obj.NLnets,obj.SortedTree_cutloop,obj.SortedCoTree_cutloop);
    [C,D,~,~,StateNamesCD]=obj.loopfixCD_num(A,B,C,D,H,s,obj.NLnets,obj.SortedTree_cutloop,obj.SortedCoTree_cutloop);

% openvar('variable_name')

% Update class
% obj.HtempAB = HtempAB;
% obj.HtempCD = HtempCD;
% obj.dependsAB = dependsAB;
% obj.savedCD = savedCD;
obj.StateNames = obj.NLnets(OrderedNamesnum,1);
obj.OutputNamesCD = StateNamesCD;
% obj.DependentNames = DependentNames;
obj.OutputNames = OutputNames;
% obj.SortedTree = SortedTree;
% obj.SortedCoTree = SortedCoTree;
obj.ConstantNames = ConstantNames;
% obj.OrderedNamesnum = OrderedNamesnum;
%obj.Codex = The_Codex;

end % That's all Folks
