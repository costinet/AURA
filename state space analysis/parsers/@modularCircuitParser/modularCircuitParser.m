classdef modularCircuitParser < NetlistCircuitParser
    %CIRCUITPARSER abstract class defining requirements of a circuitParser
    %   
    %   See also @LTspiceCircuitParser, @PLECScircuitParser
    
    properties ( SetAccess = private)

        baseParser
        moduleParser

        Nmodules = 1
    end

    % properties (SetAccess = private, Dependent)
    %     sourcefn
    %     sourcefdate
    % end

    properties (Hidden)
        sheetConnectorOutSym = '<<'
        sheetConnectorInSym = '>>'
        sheetConnectorOutFormat = @(x) startsWith(x,'<<')
        sheetConnectorInFormat = @(x) startsWith(x,'>>')
    end

    methods 
        % [Cbnd, Dbnd, hyst, switchRef] = getConstraintMatrices(obj,circuitPath)
        % loadModel(obj, fn, swseq, force)
        % updateComponentValues(obj)
        % storedTopology = saveTopology(obj,name)
        % loadTopology(obj,storedTopolpogy)


        function obj = modularCircuitParser(topology, baseFile, modularFile, Nmodules)
            obj = obj@NetlistCircuitParser(topology);

            if nargin >= 3
                obj.baseParser = NetlistCircuitParser(topology);
                obj.moduleParser = NetlistCircuitParser(topology);
    
                obj.baseParser.loadModel(baseFile);
                obj.moduleParser.loadModel(modularFile);

                obj.sourcefn = {baseFile, modularFile};
            end
            if nargin >= 4
                obj.assemble(Nmodules);
            end
        end


        function assemble(obj,Nmodules)
            obj.Nmodules = Nmodules;
            baseComps = obj.baseParser.origComponents;
            moduleComps = obj.moduleParser.origComponents;
            obj.origComponents = [];

            modNames = [{'b'} cellstr([repmat('m',Nmodules,1) num2str((1:Nmodules)')])' {'b'}];
            connNames = cat(2, cellstr(repmat('_{(',Nmodules+2,1)), ...
                join(cat(2, modNames', circshift(modNames,-1,2)' ),'-'), ...
                cellstr(repmat(')}',Nmodules+2,1)) ....
                );
            connNames = join(connNames(1:end-1,:),'');


            for i = 1:numel(baseComps)
                baseComps(i).Name = [baseComps(i).Name '(b)'];
                outNodes = obj.sheetConnectorOutFormat(baseComps(i).Nodes);
                inNodes = obj.sheetConnectorInFormat(baseComps(i).Nodes);
                zeroNodes = strcmp(baseComps(i).Nodes,'0');
                baseComps(i).Nodes(outNodes) = ...
                  cellfun(@(x) cat(2,x,connNames{1}), baseComps(i).Nodes(outNodes),'UniformOutput',false);
                baseComps(i).Nodes(inNodes) = ...
                    cellfun(@(x) cat(2,x,connNames{end}), baseComps(i).Nodes(inNodes),'UniformOutput',false);
                baseComps(i).Nodes(~outNodes & ~inNodes & ~zeroNodes) = ...
                    cellfun(@(x) cat(2,x,['_{(' modNames{1} ')}']), baseComps(i).Nodes(~outNodes & ~inNodes & ~zeroNodes),'UniformOutput',false);
                baseComps(i).Nodes = strrep(baseComps(i).Nodes,obj.sheetConnectorOutSym,'');
                baseComps(i).Nodes = strrep(baseComps(i).Nodes,obj.sheetConnectorInSym,'');
            end
            obj.origComponents = [obj.origComponents baseComps];

            for j = 1:Nmodules
                modComps = moduleComps;
                for i = 1:numel(modComps)
                    modComps(i).Name = [modComps(i).Name '(m' num2str(j) ')'];
                    outNodes = obj.sheetConnectorOutFormat(modComps(i).Nodes);
                    inNodes = obj.sheetConnectorInFormat(modComps(i).Nodes);
                    zeroNodes = strcmp(modComps(i).Nodes,'0');
                    modComps(i).Nodes(outNodes) = ...
                      cellfun(@(x) cat(2,x,connNames{j}), modComps(i).Nodes(outNodes),'UniformOutput',false);
                    modComps(i).Nodes(inNodes) = ...
                        cellfun(@(x) cat(2,x,connNames{j+1}), modComps(i).Nodes(inNodes),'UniformOutput',false);
                    modComps(i).Nodes(~outNodes & ~inNodes & ~zeroNodes) = ...
                        cellfun(@(x) cat(2,x,['_{(' modNames{j+1} ')}']), modComps(i).Nodes(~outNodes & ~inNodes & ~zeroNodes),'UniformOutput',false);
                    modComps(i).Nodes = strrep(modComps(i).Nodes,obj.sheetConnectorOutSym,'');
                    modComps(i).Nodes = strrep(modComps(i).Nodes,obj.sheetConnectorInSym,'');
                end
                obj.origComponents = [obj.origComponents modComps];
            end
            obj.sourcefn = '';
            obj.sourceType = 'modular';
            switchInds = strcmp({obj.origComponents.Type},'M') | ...
                strcmp({obj.origComponents.Type},'D');
            obj.topology.switchLabels = {obj.origComponents(switchInds).Name};
            %% Spice directives??
        end


        function storedTopology = saveTopology(obj,name)
            storedTopology = {};
            storedTopology.name = name;
            % storedTopology.components = obj.origComponents;
            % storedTopology.props = obj.netListDirectives;
            % storedTopology.switches = obj.topology.switchLabels;
            storedTopology.baseParser = saveTopology(obj.baseParser,[name '(base)']);
            storedTopology.moduleParser = saveTopology(obj.moduleParser,[name '(base)']);
        end

        function loadTopology(obj,storedTopology, Nmodules)
            obj.baseParser = NetlistCircuitParser;
            obj.moduleParser = NetlistCircuitParser;
            obj.baseParser.loadTopology(storedTopology.baseParser);
            obj.moduleParser.loadTopology(storedTopology.moduleParser);

            assemble(obj,Nmodules)
            % obj.linearizeCircuitModel2();
            obj.loadModel();

        end
    end

    methods (Hidden)
        % function [Iname, Vname] = getSwitchMeasSourceNames(obj,sName)
        % 
        % end
    end
    
end
