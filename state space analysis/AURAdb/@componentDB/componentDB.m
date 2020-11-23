classdef componentDB < handle
    %componentDB is a template for individual component databases
    
    properties (SetAccess = protected , GetAccess = protected)
        tableProps
        addlTabProps
        graphProps
        components
    end
    
    properties (Hidden, Transient, Constant)
        %SIkeys = {'f','p','n','u','m','1','k','M','G','T','P'};  
            % USE keys(SIprefixes)
        SIprefixes = containers.Map({'f','p','n',char(181),'m','','k','M','G','T'}, ...
                                    [1e-15 1e-12 1e-9 1e-6 1e-3 1e0 1e3 1e6 1e9 1e12]);        
    end
    
    properties (Hidden, Transient, Abstract, Constant)
        componentType
    end
    
    properties (Dependent)
        type
    end
    
    methods (Abstract)
        [table, addlTable, graph] = json(obj)
        sync(obj)
%         add(obj, item)
%         subsref(obj, index)
        subTable(obj, varargin) % Return a nicely-formatted table with searchable parameters
        saveDB(obj)
        loadDB(obj)
    end
    
    methods
        function obj = componentDB()
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.components = eval([class(obj.componentType) '.empty']);
        end
        
        function add(obj, item)
            assert(isa(item, class(obj.componentType)), ...
                ['class ' class(obj) ' can only be used to store objects of type ' class(obj.componentType)] );
           
            if isempty(obj.components)
                obj.components = item;
                return
            end
            
            [~, ind, ~] = intersect({obj.components.partNumber}, item.partNumber);
            if isempty(ind)
                obj.components = [obj.components item];
            else
                obj.components(ind).merge(item);
            end
        end
        
        
        function type = get.type(obj)
            type = class(obj.componentType);
        end

    end
    
    methods (Hidden)
        function l = length(obj)
            l = length(obj.components);
        end
        
        function s = size(obj)
            s = size(obj.components);
        end
        
        varargout = subsref(obj, index)
    end
        
end

