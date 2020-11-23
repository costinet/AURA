classdef transistorDB < componentDB
    %transistorDB is a database of transistor objects
    
%     Inherited properties:
%         tableProps
%         addlTabProps
%         graphProps
%         components

    properties (Dependent)
        transistors
    end
    
    properties (Hidden, Transient, Constant)
        componentType = transistor();
    end
    
    % componentDB required methods
    methods 
        [table, addlTable, graph] = json(obj)
        sync(obj)
        subTable(obj, varargin) % Return a nicely-formatted table with searchable parameters
    end
    
    methods (Hidden)

    end
    
    methods
        function obj = transistorDB()
            loadDB(obj)
        end
        
        function transistors = get.transistors(obj)
            transistors = obj.components;
        end    
        
        function saveDB(obj)
            fn = mfilename('fullpath');
            fn = [strrep(fn, '\transistorDB', '\') 'transistors.mat'];
            save(fn,'obj')
        end
        
        function loadDB(obj)
            fn = mfilename('fullpath');
            fn = [strrep(fn, '\transistorDB', '\') 'transistors.mat'];
            if isfile(fn)
                savedData = load(fn,'obj');
                obj.components = savedData.obj.transistors;
            end
        end

    end
    
end

