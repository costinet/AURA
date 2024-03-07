function components = diodeLinearSubcircuit(obj, component)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

    Vmeas = {};
    Ameas = {};
    Req = {};
    Cpar = {};
    Vf = {};

    % Note: Vmeas and Imeas are named based on what they measure, but the
    % type is the opposite.  Vm measures voltage, but is a 0A source of
    % type Im.
    Vmeas.Name = ['Vm_' component.Name]; % 0A source to measure voltage
    Ameas.Name = ['Im_' component.Name]; % 0V source to measure current
    Req.Name = ['R_' component.Name];
    Cpar.Name = ['C_' component.Name];
    Vf.Name = ['Vf_' component.Name];

    Vmeas.Type = 'Im';
    Ameas.Type = 'Vm';
    Req.Type  = 'R';
    Cpar.Type  = 'C';
    Vf.Type = 'V';

    if ~strcmp(component.Nodes{2}, '0')
        Vmeas.Nodes = component.Nodes;
        Ameas.Nodes = {component.Nodes{1}, obj.newIntNodeName(component,2,1)};
        Req.Nodes = {obj.newIntNodeName(component,2, 1), obj.newIntNodeName(component,2, 2)};
        Vf.Nodes = {obj.newIntNodeName(component,2, 2), component.Nodes{2} };
        Cpar.Nodes = component.Nodes;
    else
        Vmeas.Nodes = component.Nodes;
        Ameas.Nodes = {obj.newIntNodeName(component,1, 2), component.Nodes{2}};
        Req.Nodes = {component.Nodes{1},obj.newIntNodeName(component,1, 1)};
        Vf.Nodes = {obj.newIntNodeName(component,1, 1), obj.newIntNodeName(component,1, 2)};
        Cpar.Nodes = component.Nodes;
    end

    Vmeas.paramNames = {'V'};
    Ameas.paramNames = {'A'};
    Req.paramNames  = {'R'};
    Cpar.paramNames  = {'C'};
    Vf.paramNames  = {'Vf'};

    Vmeas.paramVals = 0;
    Ameas.paramVals = 0;
    Req.paramVals  = component.paramVals(strcmp(component.paramNames,'Roff'));
    Cpar.paramVals  = component.paramVals(strcmp(component.paramNames,'Coss') | strcmp(component.paramNames,'Cd'));
    Vf.paramVals = component.paramVals(strcmp(component.paramNames,'Vf'));

    Vmeas.paramExpressions = 0;
    Ameas.paramExpressions = 0;
    Req.paramExpressions  = component.paramExpressions(strcmp(component.paramNames,'Roff'));
    Cpar.paramExpressions  = component.paramExpressions(strcmp(component.paramNames,'Coss') | strcmp(component.paramNames,'Cd'));
    Vf.paramExpressions = component.paramExpressions(strcmp(component.paramNames,'Vf'));

    components = [Vmeas, Ameas, Req, Cpar, Vf];
end