classdef AURAdb < handle
    %AURAdb is a collection of databases used in Power Electronics Design
    
    properties (SetAccess = immutable, GetAccess = public)
        transistors
        inductors
        capacitors
        cores
        wires
        topologies
    end
    
    methods
        function obj = AURAdb()
            addpath(genpath(strrep(mfilename('fullpath'), '\AURAdb', '')));

            obj.transistors = transistorDB();
            obj.inductors = 0
            obj.capacitors = 0
            obj.cores = 0
            obj.wires = 0
        end
        
        function sync(obj)
            obj.transistors.sync();
            obj.inductors.sync();
            obj.capacitors.sync();
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

