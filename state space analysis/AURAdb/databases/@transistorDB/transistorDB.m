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
        subTable(obj, varargin) % Return a nicely-formatted table with searchable parameters
    end
    
    methods (Hidden)

    end
    
    methods
        function obj = transistorDB(skipLoad)
            if nargin == 0 || skipLoad~=1
                loadDB(obj)
            end
        end
        
        function transistors = get.transistors(obj)
            transistors = obj.components;
        end    
        
        % function saveDB(obj)
        %     fn = mfilename('fullpath');
        %     fn = [strrep(fn, '\transistorDB', '\') 'transistors.mat'];
        %     save(fn,'obj')
        % end
        % 
        % function loadDB(obj)
        %     fn = mfilename('fullpath');
        %     fn = [strrep(fn, '\transistorDB', '\') 'transistors.mat'];
        %     if isfile(fn)
        %         savedData = load(fn,'obj');
        %         obj.components = savedData.obj.transistors;
        %     end
        % end

        function clearUpdated(obj)
            for i = 1:length(obj.components)
                obj.components(i).clearUpdated();
            end
        end

        function digitizeDatasheet(obj,datasheet)
            %% 
            plotDigitizer(datasheet);
        end

    end
    
end

