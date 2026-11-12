function nl = createSimulateableNetlist(obj, simObj, options)
%CREATESIMULATEABLENETLIST Create a simulation-ready netlist.
%
% NL = CREATESIMULATEABLENETLIST(OBJ, SIMOBJ) creates a netlist
% suitable for simulation using the network information contained in
% the circuitParser OBJ and using the operation point defined by SMPSim
% object SIMOBJ
%
% NL is the resulting simulation-ready netlist, returned as text
%
% This method converts or prepares the network representation stored in the
% circuitParser for use in simulation, incorporating the operating point, 
% modulation, and adding necessary measurement sources as defined by the
% SMPSim object.
%
% See also: SMPSim, NetlistCircuitParser
arguments (Input)
    obj circuitParser
    simObj SMPSim

    options.FETKp (1,1) {mustBeNumeric} = 100000;
    options.FETVto (1,1) {mustBeNumeric} = 0.5;
    options.DiodeN (1,1) {mustBeNumeric} = 0.1;
    options.Vdr (1,1) {mustBeNumeric} = 5;
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
        modelLine = sprintf('.model %s VDMOS (Rd=%fm Rs=%fm Vto=%e Kp=%e Cgdmax=0 Cgdmin=0 Cgs=0 Cjo=%fn Rb=0 M=0)', ...
            comp.Name, comp.paramVals(1)/2*1e3, comp.paramVals(1)/2*1e3, options.FETVto, options.FETKp, comp.paramVals(end)*1e9);
        sourceLine = {['Vg' comp.Name], gateNode, comp.Nodes{2}, 'PWL REPEAT FOREVER ('};
        sourceLine = strjoin(sourceLine, ' ');

        vdsNode = makeNodeNameUnique(obj,['v_c_' comp.Name]);
        esourceLine = {['E_v_' comp.Name], vdsNode, '0', comp.Nodes{:}, '1'};
        esourceLine = [strjoin(esourceLine, ' ')];

        swvec = simObj.swvec(:,strcmp(comp.Name, simObj.switchNames));
        ts = simObj.ts;

        ted = repmat(cumsum(ts),2,1);
        %Make second one slightly bigger to have near-vertical step.
        delta = 1.5e-6*10.^floor(log10(abs(ted(2,:))));
        delta = min(delta, [diff(ted(1,:))/2, Inf]);
        ted(2,:) = ted(2,:) + delta;
        ted = [0; ted(:)];
        state =  repmat(swvec*options.Vdr,1,2)';
        state = [state(:); state(1)];

        pwl = [ted, state]';
        sourceLine = [sourceLine, sprintf('%e, %d, ', pwl(:)), ') ENDREPEAT'];

        netlist = [netlist, compLine, modelLine, sourceLine, esourceLine];
    elseif strcmp(comp.Type, 'D')
        compLine = {comp.Name, comp.Nodes{1}, comp.Nodes{2}, comp.Name};
        compLine = [strjoin(compLine, ' ')];

        Vt = 0.02585;
        Is = exp(-comp.paramVals(end)/options.DiodeN/Vt);
        modelLine = sprintf('.model %s D (N=%e Rs=%e Cjo=%e M=0 Is=%e)', ...
            comp.Name, options.DiodeN, comp.paramVals(1), comp.paramVals(3), Is);
        % modelLine = sprintf('.model %s D (Ron=%e, Roff=%e, Vfwd=%e, Vrev=%e)', ...
        %     comp.Name, comp.paramVals(1), comp.paramVals(2), comp.paramVals(end), 1e6);

        vdsNode = makeNodeNameUnique(obj,['v_c_' comp.Name]);
        esourceLine = {['E_v_' comp.Name], vdsNode, '0', comp.Nodes{:}, '1'};
        esourceLine = [strjoin(esourceLine, ' ')];

        netlist = [netlist, compLine, modelLine, esourceLine];

    elseif strcmp(comp.Type, 'C')
        name = strrep(comp.Name, '{', '');
        name = strrep(name, '}', '');
        compLine = {name, comp.Nodes{1}, comp.Nodes{2}, sprintf('%e', comp.paramVals(:))};
        compLine = [strjoin(compLine, ' ')];

        vdsNode = makeNodeNameUnique(obj,['v_' comp.Name]);
        esourceLine = {['E_v_' comp.Name], vdsNode, '0', comp.Nodes{:}, '1'};
        esourceLine = [strjoin(esourceLine, ' ')];

        netlist = [netlist, compLine, esourceLine];
    else
        name = strrep(comp.Name, '{', '');
        name = strrep(name, '}', '');
        name = [comp.Type, name(length(comp.Type)+1:end)];     %Fix for measurement sources
        compLine = {name, comp.Nodes{1}, comp.Nodes{2}, sprintf('%e', comp.paramVals(:))};
        compLine = [strjoin(compLine, ' ')];
        netlist = [netlist, compLine];
    end

end

if ~isempty(obj.netListDirectives)
    for i = 1:length(obj.netListDirectives)
        if startsWith(obj.netListDirectives(i).Name, 'K')
            dr = obj.netListDirectives(i);

            % newLines = {};
            % for j = 1:numel(dr.paramNames)
            %     line = startsWith(netlist, dr.paramNames{j});
            %     compLine = netlist{line};
            %     compNames{j} = [compLine(1:find(compLine==' ',1)-1),  '_winding'];
            %     compLine = [compNames{j}, ' ', compLine(find(compLine==' ',1):end)];
            %     match = regexp(compLine, '[+-]?(?:\d+(?:\.\d*)?|\.\d+)[eE][+-]?\d+$', 'match');
            % 
            %     if ~isempty(match)
            %         newValue = sprintf('%e', 100 * str2double(match{1}));
            %         compLine = [compLine(1:end-length(match{1})), newValue];
            %     end
            %     newLines{j} = compLine;
            % end
            % dirLine = {dr.Name, compNames{:}, dr.paramVals};

            dirLine = {dr.Name, dr.paramNames{:}, dr.paramVals};

            dirLine = [strjoin(dirLine, ' ')];
            netlist = [netlist, dirLine];%, newLines{:}];
        else
            error('netlist Directives other than K statements not implemented yet')
        end
    end
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
        line = startsWith(netlist, [simObj.stateNames{i} ' ']);
        if ~isempty(obj.netListDirectives)
            if any(any(strcmp([obj.netListDirectives(:).paramNames(:)]',simObj.stateNames{i}),2) ...
                    & strcmp([obj.netListDirectives(:).Type],'K'))
                %part of a coupled inductor, need to adjust ic to be full
                %winding current instead of only magnetizing component
    
                if any(strcmp(simObj.outputNames, ['Im_' simObj.stateNames{i}]))
                    Y = simObj.Cs(:,:,1)*simObj.Xs(:,1) + simObj.Ds(:,:,1)*simObj.u;
                    newIC = Y(strcmp(simObj.outputNames, ['Im_' simObj.stateNames{i}]));
                    netlist{line} = [netlist{line} sprintf(' ic=%e', newIC)];
                else
                    netlist{line} = [netlist{line} sprintf(' ic=%e', Xss(i,1))];
                end
            else
                netlist{line} = [netlist{line} sprintf(' ic=%e', Xss(i,1))];
            end
        else
            netlist{line} = [netlist{line} sprintf(' ic=%e', Xss(i,1))];
        end
    end
end

%%
netlist = [netlist, sprintf('.tran 0 %e 0 %e', sum(simObj.ts), min(simObj.ts)/10)];
netlist = [netlist, '.backanno'];
netlist = [netlist, '.end'];

nl = strjoin(netlist, '\n');

end