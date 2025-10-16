function [A,B,C,D,I] = solveStateSpaceRepresentation(obj)
%solveStateSpaceRepresentation solves state space matrices for
%circuitParser class.  
%   Requires that obj.components is populated.
%
%   Method adapted from L. O. Chua and P-M Lin, Compter Aided Analysis of
%   Electronic Circuits: Algorithms & Computational Techniques", Chapters 3, 6, & 8.

    assert(numel(obj.components) >= 2, 'solveStateSpaceRepresentation can only be run once a valid circuit model is added to the NetlistCircuitParser object.')

    obj.findNormalTree();
    [H,comps] = obj.hybridMatrix();
    [A,B,C,D,I,names] = obj.stateSpaceFromHybrid(H, comps);


    if isempty(obj.topology.stateLabels)
        obj.topology.stateLabels = names{1};
    else
        assert(all(strcmp(obj.topology.stateLabels,names{1})),'When parsing a new switching subinterval, it appears that the states of the circuit are different than those from a previous subinterval')
    end

    if isempty(obj.topology.outputLabels)
        obj.topology.outputLabels = names{2};
    else
        assert(all(strcmp(obj.topology.outputLabels,names{2})),'When parsing a new switching subinterval, it appears that the outputs of the circuit are different than those from a previous subinterval')
    end

    if isempty(obj.topology.inputLabels)
        obj.topology.inputLabels = names{3};
    else
        assert(all(strcmp(obj.topology.inputLabels,names{3})),'When parsing a new switching subinterval, it appears that the inputs of the circuit are different than those from a previous subinterval')
    end



end