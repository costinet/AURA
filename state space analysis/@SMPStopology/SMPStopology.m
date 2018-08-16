classdef SMPStopology < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties

        Parse
        As
        Bs
        Cs
        Ds
        
        Xi
        Bi
        
        stateLabels

        
    end
    
    methods
        function setSS(obj, As, Bs, Cs, Ds)
            obj.As = As;
            obj.Bs = Bs;
            
            if(nargin>3)
                obj.Cs = Cs;
            else
                obj.Cs = zeros(1,size(As,1));
            end
            
            if(nargin>4)
                obj.Ds = Ds;
            else
                obj.Ds = zeros(1,size(Bs,1));
            end
        end
        
        function SetABCD(obj,parse)
            obj.As = parse.Anum;
            obj.Bs = parse.Bnum;
            obj.Cs = parse.Cnum;
            obj.Ds = parse.Dnum;
        end

    end
end

