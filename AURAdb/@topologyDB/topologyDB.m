classdef topologyDB < handle
    %topologyDB is a database of topology objects
    
%     Inherited properties:
%         tableProps
%         addlTabProps
%         graphProps
%         components

    properties 
        topology
    end
    
    properties (Hidden, Transient, Constant)
        componentType = topology();
    end
    
    % componentDB required methods
    methods 
        [table, addlTable, graph] = json(obj)
        subTable(obj, varargin) % Return a nicely-formatted table with searchable parameters
    end
    
    methods (Hidden)

    end
    
    methods
        function obj = topologyDB()
            loadDB(obj)
        end
        
        
        function saveDB(obj)
            fn = mfilename('fullpath');
            fn = [strrep(fn, '\topologyDB', '\') 'topology.mat'];
            save(fn,'obj')
        end
        
        function loadDB(obj)
            fn = mfilename('fullpath');
            fn = [strrep(fn, '\topologyDB', '\') 'topology.mat'];
            if isfile(fn)
                savedData = load(fn,'obj');
                disp('Not sure what to do here ');
                
            end
        end
        
        
        function add(obj, item)

            if isempty(obj.topology)
                obj.topology = item;
                return
            end
            [~, ind, ~] = intersect({obj.topology.Name}, item.Name);
            if isempty(ind)
                obj.topology = [obj.topology item];
            else
               warning('Item not added becuase it has the same name as another topology');
            end
        end
        
        
        
    end
    
end

