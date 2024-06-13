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
        % sync(obj)
%         add(obj, item)
%         subsref(obj, index)
        subTable(obj, varargin) % Return a nicely-formatted table with searchable parameters
        % saveDB(obj)
        % loadDB(obj)
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

        function addMult(obj, items)
            assert(isa(items, class(obj.componentType)), ...
                ['class ' class(obj) ' can only be used to store objects of type ' class(obj.componentType)] );
           
            if isempty(obj.components)
                obj.components = items;
                return
            end
            
            [~, indObj, indItems] = intersect({obj.components.partNumber}, {items.partNumber});
            if isempty(indObj)
                obj.components = [obj.components items];
            else
                [~,newComps] = setdiff({items.partNumber}, {obj.components.partNumber});
                obj.components = [obj.components items(newComps)];
                for i = 1:length(indObj)
                    obj.components(indObj(i)).merge(items(indItems(i)));
                end
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

        function saveDB(obj)
            % fn = mfilename('fullpath');
            % fn = [strrep(fn, '\transistorDB', '\') 'transistors.mat'];
            % save(fn,'obj')
            DBpath = fullfile(userpath, 'SMPSToolbox', 'AURAdb', [class(obj.componentType) 'DB'],filesep);
            if isempty(dir(DBpath))
                mkdir(DBpath)
            end
            save(fullfile(DBpath,[class(obj.componentType) 's.mat']),'obj');
        end
        
        function loadDB(obj, reload)
            if nargin == 1
                reload = 0;
            end

            DBpath = fullfile(userpath, 'SMPSToolbox', 'AURAdb', [class(obj.componentType) 'DB'],filesep);
            if isempty(dir(DBpath))
                mkdir(DBpath)
                return
            end
            if ~isempty(dir(DBpath))
                matfiles = dir(fullfile(DBpath, '*.mat'));
                if ~isempty(matfiles)
                    for i = 1:numel(matfiles)
                        savedData = load(fullfile(DBpath,matfiles(i).name),'obj');
                        % for j = 1:numel(savedData.obj.components)
                        %     obj.add(savedData.obj.components(j));
                        % end
                        obj.addMult(savedData.obj.components);
                    end
                    if numel(matfiles) > 1
                        warning(['Only one saved database should be present in ' DBpath '.  Possible duplication of data may occur']);
                    end
                end
                if isempty(matfiles) || reload
                    libfiles = dir(fullfile(DBpath,'lib','*.mat'));
                    if ~isempty(libfiles)
                        for i = 1:numel(libfiles)
                            savedData = load(fullfile(DBpath,'lib',libfiles(i).name),'obj');
                            % for j = 1:numel(savedData.obj.components)
                            %     obj.add(savedData.obj.components(j));
                            % end
                            if isa(savedData.obj, 'componentDB')
                                obj.addMult(savedData.obj.components);
                            elseif isa(savedData.obj, 'component')
                                obj.addMult(savedData.obj);
                            else
                                error(['Incorrect formatting of file ' fullfile(DBpath,'lib',libfiles(i).name)]);
                            end
                        end
                    end
                end
            else
                warning(['Unable to access user documents path at ' userpath]);
            end
        end
        
        function ind = end(obj, k, n)
            %% redefined to aid in subsref
            ind = length(obj.components);
        end
        
        function n = numArgumentsFromSubscript(obj,s,indexingContext)
            %% redefined to aid in subsref
            if strcmp(s(1).type, '()') 
                if strcmp(s(1).subs{:}, ':')
                    n = length(obj.components(s(1).subs{:}));
                else
                    if length(s) == 1
                        n = builtin('numArgumentsFromSubscript',obj.components,s,indexingContext);
                    else
                        n = 0;
                        comps = obj.components(s(1).subs{:});
                        for i = 1:length(comps)
                            n = n + builtin('numArgumentsFromSubscript',obj.components(i),s(2:end),indexingContext);
                        end
                    end
                end
            else
                n = builtin('numArgumentsFromSubscript',obj,s,indexingContext);
            end
        end
        
%         varargout = subsref(obj, index)
    end
        
end

