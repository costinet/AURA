function components = XFdependentSourceSubcircuit(obj, directive, components)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

    
    for i = 1:length(directive.paramNames)
        coupledIndLocs(i) = find(strcmp({components.Name}, directive.paramNames{i}));
        coupledInductors(i) = components(coupledIndLocs(i));
    end
    
    k = eval(directive.paramVals);

    if k==1
        %%Use ideal transformer model
        Vi = 2; % for the moment, assume the first inductor is the voltage reference.  
                % later, should look at impedance seen by each coil
        N(Vi) = 1;
        LV = coupledInductors(Vi).paramVals(strcmp(coupledInductors(Vi).paramNames,'L'));
        isource = {};
        isource.Name = ['F_' coupledInductors(Vi).Name];
        isource.Type = 'F';
        isource.Nodes = coupledInductors(Vi).Nodes;
        isource.paramNames = {};
        isource.paramVals = {};


        for i = setdiff(1:length(coupledInductors),Vi)
            Li = coupledInductors(i).paramVals(strcmp(coupledInductors(i).paramNames,'L'));
            N(i) = sqrt(Li/LV);
            
            Vsource = {};
            Vsource.Name = ['E_' coupledInductors(i).Name];
            Vsource.Type = 'E';
            Vsource.Nodes = [coupledInductors(i).Nodes];
            Vsource.paramNames = {'Inam', 'gain'};   %This is not valid LTspice, but it is easier to make it a two-node component
            Vsource.paramVals = {isource.Name, N(i)};

            isource.paramNames = [isource.paramNames, ...
                {['Vnam' num2str(i)], ['gain' num2str(i)]}];
            isource.paramVals = [isource.paramVals, ...
                {Vsource.Name, N(i)}];

            components(coupledIndLocs(i)) = Vsource;    %Replace all secondaries

        end


        components = [components(1:Vi) isource components(Vi+1:end)];   %Leave primary (magnetizing) inductance in
    else
        error('coupled inductors with k=/=1 not currently supported');
    end

end