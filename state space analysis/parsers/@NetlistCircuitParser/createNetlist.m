function nl = createNetlist(obj, simObj)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    obj circuitParser
    simObj SMPSim
end

netlist = {'*Test Netlist'};

for i = 1:numel(obj.origComponents)
    comp = obj.origComponents(i);
    if strcmp(comp.Type, 'M')
        % Add component to netlist
            %Here we add a new gate node, named for the FET
        gateNode = makeNodeNameUnique(obj,['g' comp.Name]);
        compLine = {comp.Name, comp.Nodes{1}, [gateNode], comp.Nodes{2}, comp.Nodes{2}, comp.Name};
        compLine = [strjoin(compLine, ' ')];
        modelLine = sprintf('.model %s VDMOS (Rd=%fm Rs=%fm Vto=0.5 Kp=100k Cgdmax=0 Cgdmin=0 Cgs=0 Cjo=%fn Rb=0)', ...
            comp.Name, comp.paramVals(1)/2*1e3, comp.paramVals(1)/2*1e3, comp.paramVals(end)*1e9);
        sourceLine = {['Vg' comp.Name], gateNode, comp.Nodes{2}, 'PWL REPEAT FOREVER ('};
        sourceLine = strjoin(sourceLine, ' ');

        vdsNode = makeNodeNameUnique(obj,['v_c_' comp.Name]);
        esourceLine = {['E_v_' comp.Name], vdsNode, '0', comp.Nodes{:}, '1'};
        esourceLine = [strjoin(esourceLine, ' ')];

        swvec = simObj.swvec(:,strcmp(comp.Name, simObj.switchNames));
        ts = simObj.ts;

        ted = repmat(cumsum(ts),2,1);
        %Make second one slightly bigger to have near-vertical step.
        delta = 1e-6*10.^floor(log10(abs(ted(2,:))));
        ted(2,:) = ted(2,:) + delta;
        ted = [0; ted(:)];
        state =  repmat(swvec*5,1,2)';
        state = [state(:); state(end)];

        pwl = [ted, state]';
        sourceLine = [sourceLine, sprintf('%e, %d, ', pwl(:)), ') ENDREPEAT'];

        netlist = [netlist, compLine, modelLine, sourceLine, esourceLine];
    elseif strcmp(comp.Type, 'D')
        compLine = {comp.Name, comp.Nodes{1}, comp.Nodes{2}, comp.Name};
        compLine = [strjoin(compLine, ' ')];
        modelLine = sprintf('.model %s D (N=.001, Rs=%e, Cjo=%e)', ...
            comp.Name, comp.paramVals(1), comp.paramVals(end));

        netlist = [netlist, compLine, modelLine];
    else
        compLine = {comp.Name, comp.Nodes{1}, comp.Nodes{2}, sprintf('%e', comp.paramVals(:))};
        compLine = [strjoin(compLine, ' ')];
        netlist = [netlist, compLine];
    end

end

if ~isempty(obj.netListDirectives)
    error('netlist Directives not implemented yet')
end

%% Initial conditions
Xss = simObj.Xs;
if isempty(Xss)
    Xss = simObj.steadyState;
end
for i = 1:numel(simObj.stateNames)
    if startsWith(simObj.stateNames(i), 'C_M') || startsWith(simObj.stateNames(i), 'C_D')
        fet = obj.origComponents(strcmp({obj.origComponents.Name}, simObj.stateNames{i}(3:end)));
        compLine = {simObj.stateNames{i}, fet.Nodes{:}, '0', sprintf('ic=%e', Xss(i,1))};
        compLine = [strjoin(compLine, ' ')];
        netlist = [netlist, compLine];
    elseif startsWith(simObj.stateNames(i), 'L') || startsWith(simObj.stateNames(i), 'C')
        line = startsWith(netlist, simObj.stateNames(i));
        netlist{line} = [netlist{line} sprintf(' ic=%e', Xss(i,1))];
    end
end

%%
netlist = [netlist, sprintf('.tran %e', sum(simObj.ts))];
netlist = [netlist, '.backanno'];
netlist = [netlist, '.end'];

nl = strjoin(netlist, '\n');

end