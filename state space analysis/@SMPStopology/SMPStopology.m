classdef SMPStopology < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        As
        Bs
        Cs
        Ds
        
        Xi
        Bi
        
        stateLabels
        outputLabels
        inputLabels

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
    end
    
end

