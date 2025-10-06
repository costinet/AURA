classdef (Abstract) componentDB < handle
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

    properties (Hidden)
        modified = 0
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
            %componentDB Construct an instance of this class
            %   componentDB is a generic class that is a template for
            %   components of a specific type (e.g. transistor, capacitor,
            %   etc.)
            obj.components = eval([class(obj.componentType) '.empty']);
        end
        
        function add(obj, item)
            assert(isa(item, class(obj.componentType)) || isa(item, obj.namespaceFreeComponentType), ...
                ['class ' class(obj) ' can only be used to store objects of type ' class(obj.componentType)] );
           
            if isempty(obj.components)
                if isa(item, class(obj.componentType))
                    obj.components = item;
                    obj.modified = 1;
                else

                end
                return
            end
            
            if isempty({obj.components.partNumber}) || isempty(item.partNumber)
                ind = [];
            else
                [~, ind, ~] = intersect({obj.components.partNumber}, item.partNumber);
            end
            if isempty(ind)
                obj.components = [obj.components item];
                obj. modified  = 1;
            else
                obj.components(ind).merge(item);
                if obj.components(ind).upDated == 1
                    obj.modified = 1;
                end

            end
        end

        function addMult(obj, items)
            assert(isa(items, class(obj.componentType)) || isa(items, obj.namespaceFreeComponentType), ...
                ['class ' class(obj) ' can only be used to store objects of type ' class(obj.componentType)] );
           
            if isempty(obj.components)
                if strcmp(class(items), class(obj.componentType))
                    %using strcmp instead of isa to ignore superclass
                    obj.components = items;
                    obj.modified = 1;
                else
                    error('This should never be reached')
                    % This block is reached when we need to convert from
                    % pre-namespace databases
                    % This should no longer be neceesary with the shim
                    % classes in back-compat
                    % pathFolders = strsplit(path, pathsep);
                    % rfToolboxFolders = [pathFolders(contains(pathFolders, '\toolbox\rf\', 'IgnoreCase', true))...
                    %     pathFolders(contains(pathFolders, '\toolbox\rfpcb\', 'IgnoreCase', true))];
                    % if ~isempty(rfToolboxFolders)
                    %      rmpath(rfToolboxFolders{:});
                    % end
                    % for i = 1:length(items)
                    %     comp = feval(class(obj.componentType));
                    %     if isa(comp,'smps.component')
                    %         comp.partNumber = items(i).partNumber;
                    %         merge(comp, items(i));
                    %     elseif isa(comp,'smps.components.storedTopology')
                    %         obj.add(items(i));
                    %     end
                    % 
                    %     obj.components = [obj.components comp];
                    % end
                    %  if ~isempty(rfToolboxFolders)
                    %      addpath(rfToolboxFolders{:});
                    % end
                end
                return
            end
            
            [~, indObj, indItems] = intersect({obj.components.partNumber}, {items.partNumber});
            if isempty(indObj)
                obj.components = [obj.components items];
                obj.modified = 1;
            else
                [~,newComps] = setdiff({items.partNumber}, {obj.components.partNumber});
                obj.components = [obj.components items(newComps)];
                for i = 1:length(indObj)
                    obj.components(indObj(i)).merge(items(indItems(i)));
                    if obj.components(indObj(i)).upDated == 1
                        obj.modified = 1;
                    end
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

        % Moved to AURAdb
        % function delete(obj)
        %     % delete is overwritted for componentDB classes to check if
        %     % anything has been modified and give the user a warning to
        %     % prevent data loss.
        %     % 
        %     if obj.modified == 1
        %         selection = questdlg(['You are about to delete an object of type ' class(obj) ' which has been modified since loading. Would you like to save the database before deleting?'], ...
        %             "Confirm Deleting Unsaved Data");
        %         % outfn = websave(fullfile(userpath, 'SMPSToolbox', 'AURAdb', [latestVersion '.zip']) , releaseURL);
        %         if strcmp(selection, 'Yes')
        %             obj.saveDB();
        %         end
        %     end
        % end

        function saveDB(obj, alreadyConfirmed)
            % saveDB saves database to mat file in user's folder
            %   

            if nargin == 1
                alreadyConfirmed = 0;
            end

            % fn = mfilename('fullpath');
            % fn = [strrep(fn, '\transistorDB', '\') 'transistors.mat'];
            % save(fn,'obj')
            % DBpath = fullfile(userpath, 'SMPSToolbox', 'AURAdb', [obj.namespaceFreeComponentType 'DB'],filesep);
            DBpath = fullfile(userpath, 'SMPSToolbox', 'AURAdb', [class(obj.componentType) 'DB'],filesep);
            if isempty(dir(DBpath))
                mkdir(DBpath)
                alreadyConfirmed = 1;
            end
            
            
            fn = fullfile(DBpath,[class(obj.componentType) 's.mat']);
            if exist(fn,'file') && ~alreadyConfirmed
                selection = questdlg(['Are you sure you want to overwrite the ' class(obj.componentType) ' database?'], ...
                    "Confirm Overwriting Database");
                if strcmp(selection, 'Yes')
                    alreadyConfirmed = 1;
                else
                    alreadyConfirmed = 0;
                    return
                end
            end

            if alreadyConfirmed == 1
                for i = 1:length(obj.components)
                    obj.components(i).resetUpDated();
                    obj.modified = 0;
                end
                save(fn,'obj');
            end
            
            
        end
        
        function loadDB(obj, reload)
            if nargin == 1
                reload = 0;
            end

            DBpath = fullfile(userpath, 'SMPSToolbox', 'AURAdb', [class(obj.componentType) 'DB'],filesep);
            if isempty(dir(DBpath))
                % mkdir(DBpath)
                % Check if old version exists
                DBpath  = fullfile(userpath, 'SMPSToolbox', 'AURAdb', [obj.namespaceFreeComponentType 'DB'],filesep);
                if isempty(dir(DBpath))
                    return
                else
                    loadOutdatedDatabase(obj, DBpath)
                    obj.modified = 1;
                    return
                end
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
                            if isa(savedData.obj, 'smps.componentDB')
                                obj.addMult(savedData.obj.components);
                            elseif  isa(savedData.obj, 'smps.component')
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
           obj.modified = 0;
        end

        function loadOutdatedDatabase(obj, DBpath)
            % Back-compatability function for loading databases that were
            % saved without the smps namespace.  
            %   loadOutdatedDatabase follows loadDB except that it will
            %   first try to remove the rf toolbox, if present, due to
            %   naming conflicts, then allow a broader range of classes to
            %   be loaded. The function restores the rf toolbox after
            %   completion or in the event of an error.

            pathFolders = strsplit(path, pathsep);
            rfToolboxFolders = [pathFolders(contains(pathFolders, '\toolbox\rf\', 'IgnoreCase', true))...
                pathFolders(contains(pathFolders, '\toolbox\rfpcb\', 'IgnoreCase', true))];
            if ~isempty(rfToolboxFolders)
                 rmpath(rfToolboxFolders{:});
            end

            try
                if ~isempty(dir(DBpath))
                    matfiles = dir(fullfile(DBpath, '*.mat'));
                    if ~isempty(matfiles)
                        for i = 1:numel(matfiles)
                            savedData = load(fullfile(DBpath,matfiles(i).name),'obj');
                            obj.addMult(savedData.obj.components);
                        end
                        if numel(matfiles) > 1
                            warning(['Only one saved database should be present in ' DBpath '.  Possible duplication of data may occur']);
                        end
                    end
                    if isempty(matfiles) 
                        libfiles = dir(fullfile(DBpath,'lib','*.mat'));
                        if ~isempty(libfiles)
                            for i = 1:numel(libfiles)
                                savedData = load(fullfile(DBpath,'lib',libfiles(i).name),'obj');
                                % for j = 1:numel(savedData.obj.components)
                                %     obj.add(savedData.obj.components(j));
                                % end
                                if isa(savedData.obj, 'componentDB') || isa(savedData.obj, 'capacitorDB') || isa(savedData.obj, 'transistorDB') ||  isa(savedData.obj, 'smps.componentDB')
                                    obj.addMult(savedData.obj.components);
                                elseif isa(savedData.obj, 'component') || isa(savedData.obj, 'capacitor') || isa(savedData.obj, 'transistor') ||  isa(savedData.obj, 'smps.component')
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

                addpath(rfToolboxFolders{:});
            catch e
                addpath(rfToolboxFolders{:});
                rethrow(e);
            end

        end
        
        function ind = end(obj, k, n)
            %% redefined to aid in subsref
            ind = length(obj.components);
        end

        function tf = isempty(obj)
            % Overload to pass to components array
            tf = isempty(obj.components);
        end

        function compName = namespaceFreeComponentType(obj)
             compName = split(class(obj.componentType),'.');
             compName = compName{end};
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

