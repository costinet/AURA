function components = switchLinearSubcircuit(obj, component)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

    Vmeas = {};
    Ameas = {};
    Req = {};
    Cpar = {};

    % Note: Vmeas and Imeas are named based on what they measure, but the
    % type is the opposite.  Vm measures voltage, but is a 0A source of
    % type Im.
    Vmeas.Name = ['Vm_' component.Name]; % 0A source to measure voltage
    Ameas.Name = ['Im_' component.Name]; % 0V source to measure current
    Req.Name = ['R_' component.Name];
    Cpar.Name = ['C_' component.Name];

    Vmeas.Type = 'Im';
    Ameas.Type = 'Vm';
    Req.Type  = 'R';
    Cpar.Type  = 'C';

    if ~strcmp(component.Nodes{2}, '0')
        Vmeas.Nodes = component.Nodes;
        Ameas.Nodes = {component.Nodes{1}, obj.newIntNodeName(component,2, 1)};
        Req.Nodes = {obj.newIntNodeName(component,2, 1), component.Nodes{2}};
        Cpar.Nodes = component.Nodes;
    else
        Vmeas.Nodes = component.Nodes;
        Ameas.Nodes = {obj.newIntNodeName(component,1, 1), component.Nodes{2}};
        Req.Nodes = {component.Nodes{1}, obj.newIntNodeName(component,1, 1)};
        Cpar.Nodes = component.Nodes;
    end

    Vmeas.paramNames = {'V'};
    Ameas.paramNames = {'A'};
    Req.paramNames  = {'R'};
    Cpar.paramNames  = {'C'};

    Vmeas.paramVals = 0;
    Ameas.paramVals = 0;
    Req.paramVals  = component.paramVals(strcmp(component.paramNames,'Roff'));
    Cpar.paramVals  = component.paramVals(strcmp(component.paramNames,'Coss') | strcmp(component.paramNames,'Cd'));

    Vmeas.paramExpressions = {};
    Ameas.paramExpressions = {};
    Req.paramExpressions  = component.paramExpressions(strcmp(component.paramNames,'Roff'));
    Cpar.paramExpressions  = component.paramExpressions(strcmp(component.paramNames,'Coss') | strcmp(component.paramNames,'Cd'));

    %% Deal with zero (neglected) parasitics
    components = [Vmeas, Ameas, Req]; 
    if Cpar.paramVals > 0  
        components = [components, Cpar];
    end
end