classdef AURAdb < handle
    %% AURAdb is a collection of databases used in Power Electronics Design
    %
    % See also SMPSim, @transistor, @capacitor
    
    properties (SetAccess = immutable, GetAccess = public)
        transistors
        inductors
        capacitors
        cores
        wires
        topologies
    end
    
    methods
        function obj = AURAdb(softLoad)

            % addpath(genpath(strrep(mfilename('fullpath'), '\AURAdb', '')));
            if nargin == 1 && softLoad
                return
            end

            obj.transistors = transistorDB();
            obj.inductors = 0;
            obj.capacitors = capacitorDB();
            obj.cores = 0;
            obj.wires = 0;
        end
        
        function sync(obj)
            obj.transistors.sync();
            obj.inductors.sync();
            obj.capacitors.sync();
        end

        function updateLibraries(obj)
            gitRelease = webread('https://api.github.com/repos/costinet/AURAdb/releases/latest');
            latestVersion = gitRelease.tag_name;
            releaseURL = gitRelease.zipball_url;
            if ~exist(fullfile(userpath, 'SMPSToolbox', 'AURAdb', latestVersion), 'dir')
                disp(['Downloading latest libraries from ' releaseURL]);
                % outfn = websave(fullfile(userpath, 'SMPSToolbox', 'AURAdb', [latestVersion '.zip']) , releaseURL);
                fns = unzip(releaseURL, fullfile(userpath, 'SMPSToolbox', 'AURAdb', latestVersion));
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
            disp(['Cleaning up']);
            if exist(fullfile(userpath, 'SMPSToolbox', 'AURAdb', latestVersion), 'dir')
                rmdir(fullfile(userpath, 'SMPSToolbox', 'AURAdb', latestVersion),'s');
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

