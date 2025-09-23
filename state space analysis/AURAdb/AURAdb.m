classdef AURAdb < handle
    %% AURAdb is a collection of databases used in Power Electronics Design
    %
    % See also SMPSim, @transistor, @capacitor
    
    properties (SetAccess = immutable, GetAccess = public)
        transistors = smps.databases.transistorDB(1)
        capacitors = smps.databases.capacitorDB(1)
        % inductors
        % cores
        % wires
        topologies = smps.databases.topologyDB(1)

    end

    properties (Dependent, Hidden)
        allComponents
    end
    
    methods
        function obj = AURAdb(softLoad)

            % addpath(genpath(strrep(mfilename('fullpath'), '\AURAdb', '')));
            if nargin == 1 && softLoad
                return
            end

            obj.transistors = smps.databases.transistorDB();
            obj.capacitors = smps.databases.capacitorDB();
            obj.topologies = smps.databases.topologyDB();
            
            % try
                % obj.capacitors = smps.databases.capacitorDB();
            % catch e
            %     % RF Toolbox also has a capacitor() function.  Check if it
            %     % is conflicting and report to user.
            %     if endsWith(e.stack(1).file, fullfile('rf','capacitor.m'))
            %         e = e.addCause(MException('CAPACITOR:shadowedByRFToolbox',[ ...
            %             'function capacitor() for AURAdb components is shadowed by a identically-named function in the RF Toolbox. ' ...
            %             'Move the AURAdb\\components\\capacitor folder lower than the RF toolbox in the matlab path']));
            %         throw(e);
            %     else
            %         rethrow(e)
            %     end
            % 
            % end
            % obj.inductors = 0;
            % obj.cores = 0;
            % obj.wires = 0;

            if isempty(obj.transistors) && isempty(obj.capacitors) && isempty(obj.topologies)
                warning('AURAdb is empty.  Run AURAdb(1).updateLibraries() to sync the lates libraries from the repository')
            end
        end

        function add(obj, component)
            if isscalar(component)
                for i = 1:length(obj.allComponents)
                    if isa(component, class(obj.allComponents{i}.componentType))
                        add(obj.allComponents{i}, component);
                        return
                    end
                end
            else
                for i = 1:length(obj.allComponents)
                    if isa(component(1), class(obj.allComponents{i}.componentType))
                        obj.allComponents{i}.addMult(component);
                        return
                    end
                end
            end
            e = MException('AURAdb:InvalidComponent', 'Invalid object being added to AURAdb.  Inputs to the add() function must be components of the same class as one of the component datbases.');
            throw(e);

        end
        
        % function sync(obj)
        %     obj.transistors.sync();
        %     obj.inductors.sync();
        %     obj.capacitors.sync();
        % end

        function updateLibraries(obj)
            gitRelease = webread('https://api.github.com/repos/costinet/AURAdb/releases/latest');
            latestVersion = gitRelease.tag_name;
            releaseURL = gitRelease.zipball_url;
            if ~exist(fullfile(userpath, 'SMPSToolbox', 'AURAdb', latestVersion), 'dir')
                disp(['Downloading latest libraries from ' releaseURL]);
                selection = questdlg(['AURAdb device library version ' latestVersion ' is available from GitHub.  Would you like to install?'], ...
                    "Confirm Install");
                % outfn = websave(fullfile(userpath, 'SMPSToolbox', 'AURAdb', [latestVersion '.zip']) , releaseURL);
                if strcmp(selection, 'Yes')
                    fns = unzip(releaseURL, fullfile(userpath, 'SMPSToolbox', 'AURAdb', latestVersion));
                else
                    fns = [];
                    warning('Run updateToolbox() anytime to update both AURA and AURAdb.');
                end
            else
                disp('AURAdb is updated to latest version');
                fns = [];
            end
            if ~isempty(fns)
                disp('Extracting files');
                for db = {'capacitorDB', 'transistorDB', 'inductorDB'}
                    dirLoc = endsWith(fns,[db{1} filesep]);
                    if any(dirLoc)
                       dirLoc = find(dirLoc,1);
                       movefile(fns{dirLoc}, fullfile(userpath, 'SMPSToolbox', 'AURAdb', filesep));
                    end
                end
            end
            disp('Cleaning up');
            if exist(fullfile(userpath, 'SMPSToolbox', 'AURAdb', latestVersion), 'dir')
                % rmdir(fullfile(userpath, 'SMPSToolbox', 'AURAdb', latestVersion),'s');
                % mkdir(fullfile(userpath, 'SMPSToolbox', 'AURAdb', latestVersion));
            end
            disp('AURAdb Sync completed');
        end

        function comps = get.allComponents(obj)
            comps = {obj.transistors, obj.capacitors, obj.topologies};
       end
    end

   methods (Hidden)

       function delete(obj)
            % delete is overwritted for componentDB classes to check if
            % anything has been modified and give the user a warning to
            % prevent data loss.
            % 
            modified = cellfun(@(x) x.modified,obj.allComponents);
            if any(modified)
                selection = questdlg('AURAdb contains modified databases. Would you like to save the database before deleting?', ...
                        "Confirm Deleting Unsaved Data");
            else
                return
            end
            if strcmp(selection, 'Yes')
                for db = obj.allComponents
                    if db.modified == 1
                         db.saveDB(1);
                    end
                end
            end
        end
    end
    
    methods (Static, Hidden)
        function hash = SHA256(filename)
            text = fileread(filename);
            md = java.security.MessageDigest.getInstance('SHA-256');
            hash = sprintf('%2.2x', typecast(md.digest(uint8(text)), 'uint8')');
        end
    end


end

