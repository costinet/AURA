classdef transistorDB < componentDB
    %transistorDB Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %% Inherited
%         componentType
%         relDB
%         nrDB
%         relDBfields
    end
    
    methods
        function obj = transistorDB(relfn,nrfn)
            %transistorDB Construct an instance of this class
            %   Detailed explanation goes here
            if ~isfile(relfn)
                conn = sqlite(relfn,'create');
                relDBinit = "create table FETs (partno varchar(255) NOT NULL UNIQUE, ron NUMERIC, vds NUMERIC, Id NUMERIC, Qg NUMERIC, Rjc Numeric)";
                exec(conn, relDBinit)
                close(conn);
            end
            if ~isfile(nrfn)
                fid = fopen(nrfn,'w');
                fwrite(fid, '[]');
                fclose(fid);
            end
            obj@componentDB(relfn,nrfn); % Call superclass constructor                   
            obj.componentType = 'Transistor';
            
            getcols = ["SELECT * FROM FETs "];
            [results,metadata] = select(obj.relDB, getcols);
            
            obj.relDBfields = metadata.Properties
            sqlread(obj.relDB,'FETs')
        end
        
        function addTransistor(obj,transistor)
            %addTransistor Summary of this method goes here
            %   Detailed explanation goes here
            insertFET = ['INSERT into FETs (partno, ron, vds, Id, Qg, Rjc) VALUES ("', ...
                transistor.partno, '", ' ...
                num2str(transistor.ron,'%f'), ', ' ...
                num2str(transistor.Vbr,'%f'), ', ' ...
                num2str(transistor.Id,'%f'), ', ' ...
                num2str(transistor.qg,'%0.15f'), ', ' ...
                num2str(transistor.Rjc,'%f'), ') ' ...
                ];
            exec(obj.relDB, insertFET);
        end
    end
end

