classdef topologyDB < componentDB
    %topologyDB is a database of stored SMPS converter topology objects
    
%     Inherited properties:
%         tableProps
%         addlTabProps
%         graphProps
%         components

    properties (Dependent)
        topologies
    end
    
    properties (Hidden, Transient, Constant)
        componentType = storedTopology();
    end
    
    % componentDB required methods
    methods 
        [table, addlTable, graph] = json(obj)
        subTable(obj, varargin) % Return a nicely-formatted table with searchable parameters
    end
    
    
    methods
        function obj = topologyDB(skipLoad)
            if nargin == 0 || skipLoad~=1
                loadDB(obj)
            end
        end
        
        function topologies = get.topologies(obj)
            topologies = obj.components;
        end    

        function clearUpdated(obj)
            for i = 1:length(obj.components)
                obj.components(i).clearUpdated();
            end
        end

        function add(obj, item, overwrite)
            arguments
                obj topologyDB
                item storedTopology
                overwrite = 0
            end
            assert(isa(item, class(obj.componentType)), ...
                ['class ' class(obj) ' can only be used to store objects of type ' class(obj.componentType)] );
           
            if isempty(obj.components)
                obj.components = item;
                return
            end
            
            [~, ind, ~] = intersect({obj.components.name}, item.name);
            if isempty(ind)
                obj.components = [obj.components item];
            else
                if overwrite
                    obj.components(ind) = item;
                else
                    warning(['Topology with name ' item.name ' already exists in the topologyDB']);
                end
            end
        end

 

    end

    methods (Hidden)
        % function remove(obj,ind)
        %     obj.components(ind) = [];
        % end
    end
    
end

