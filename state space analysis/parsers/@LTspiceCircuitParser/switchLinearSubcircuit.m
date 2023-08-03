function components = switchLinearSubcircuit(obj, component)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

    Vmeas = {};
    Ameas = {};
    Req = {};
    Cpar = {};

    Vmeas.Name = ['Im_' component.Name];
    Ameas.Name = ['Vm_' component.Name];
    Req.Name = ['R_' component.Name];
    Cpar.Name = ['C_' component.Name];

    Vmeas.Type = 'V';
    Ameas.Type = 'I';
    Req.Type  = 'R';
    Cpar.Type  = 'C';

    if ~strcmp(component.Nodes{2}, '0')
        Vmeas.Nodes = component.Nodes;
        Ameas.Nodes = {component.Nodes{1}, [component.Nodes{2} 'm']};
        Req.Nodes = {[component.Nodes{2} 'm'], component.Nodes{2}};
        Cpar.Nodes = component.Nodes;
    else
        Vmeas.Nodes = component.Nodes;
        Ameas.Nodes = {[component.Nodes{1} 'm'], component.Nodes{2}};
        Req.Nodes = {component.Nodes{1}, [component.Nodes{1} 'm']};
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

    components = [Vmeas, Ameas, Req, Cpar];
end