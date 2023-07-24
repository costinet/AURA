classdef componentDB < handle
    %componentDB Summary of this class goes here
    %   Detailed explanation goes here
    
    properties 
        componentType
        relDB
        nrDB
        relDBfields
    end
    
    %% Private Methods
    methods (Access = private)
        function openRelDB(obj, relfn)
            %openRelDB load relational database
            %   Loads all data from json-encoded non-relational database
            %   and decodes data into a MATLAB struct
            fid = fopen(relfn);
            obj.relDB = fid;
            
%             [filepath,name,ext] = fileparts(filename)
            
%            conn = sqlite(relfn); 
%            obj.relDB = conn;
        end
        
        function openNRDB(obj, nrfn)
            %openNRDB load non-relational database
            %   Loads all data from json-encoded non-relational database
            %   and decodes data into a MATLAB struct
            edit nrfn;
%            obj.nrDB = jsondecode(fileread(nrfn));
        end
  
    end
    
    methods
        %% Constructors
        function obj = componentDB(relfn,nrfn)
            %componentDB Construct an instance of this class
            %   obj = componentDB(relfn,nrfn)
            openRelDB(obj,relfn);
            openNRDB(obj,nrfn);

        end
        
        %% Methods
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
        
        %% Destructors
        function delete(obj)
           close(obj.relDB);
        end
    end
end

