classdef AURAdb < handle
    %% AURAdb is a collection of databases used in Power Electronics Design
    %
    % See also SMPSim, @transistor, @capacitor
    
    properties (SetAccess = immutable, GetAccess = public)
        transistors
        capacitors
        % inductors
        % cores
        % wires
        topologies
    end
    
    methods
        function obj = AURAdb(softLoad)

            % addpath(genpath(strrep(mfilename('fullpath'), '\AURAdb', '')));
            if nargin == 1 && softLoad
                return
            end

            obj.transistors = transistorDB();
            
            try
                obj.capacitors = capacitorDB();
            catch e
                % RF Toolbox also has a capacitor() function.  Check if it
                % is conflicting and report to user.
                if endsWith(e.stack(1).file, fullfile('rf','capacitor.m'))
                    e = e.addCause(MException('CAPACITOR:shadowedByRFToolbox',[ ...
                        'function capacitor() for AURAdb components is shadowed by a identiclaly-named function in the RF Toolbox. ' ...
                        'Move the AURAdb\\components\\capacitor folder lower than the RF toolbox in the matlab path']));
                    throw(e);
                else
                    rethrow(e)
                end

            end
            % obj.inductors = 0;
            % obj.cores = 0;
            % obj.wires = 0;

            obj.topologies = topologyDB();
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
                % outfn = websave(fullfile(userpath, 'SMPSToolbox', 'AURAdb', [latestVersion '.zip']) , releaseURL);
                fns = unzip(releaseURL, fullfile(userpath, 'SMPSToolbox', 'AURAdb', latestVersion));
            else
                disp('AURAdb already updated to lates version');
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

