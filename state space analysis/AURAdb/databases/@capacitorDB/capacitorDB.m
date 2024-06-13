classdef capacitorDB < componentDB
    %capacitorDB is a database of capacitor objects
    
%     Inherited properties:
%         tableProps
%         addlTabProps
%         graphProps
%         components

    properties (Dependent)
        capacitors
    end
    
    properties (Hidden, Transient, Constant)
        componentType = capacitor();
    end
    
    % componentDB required methods
    methods 
        [table, addlTable, graph] = json(obj)
        subTable(obj, varargin) % Return a nicely-formatted table with searchable parameters
    end
    
    methods (Hidden)

    end
    
    methods
        function obj = capacitorDB(skipLoad)
            if nargin == 0 || skipLoad~=1
                loadDB(obj)
            end
        end
        
        function capacitors = get.capacitors(obj)
            capacitors = obj.components;
        end    
        

        function clearUpdated(obj)
            for i = 1:length(obj.components)
                obj.components(i).clearUpdated();
            end
        end



    end
    
end

