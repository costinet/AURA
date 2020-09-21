classdef SMPStopology < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties

        Parser
        As
        Bs
        Cs
        Ds
        
        Asym
        Bsym
        Csym
        Dsym
        
        filename % The filename that contains the topology for this covnerter
        
        stateLabels={}
        stateLabels_Opp={}
        Xi
        Bi
        outputLabels
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

        function []=parse(obj,filename)
            % confrim that pasre class is empty??!!??
            obj.filename = filename;
            obj.Parser.initialize(filename);
            obj.Parser.ABCD;
            obj.Asym = obj.Parser.Asym;
            obj.Bsym = obj.Parser.Bsym;
            obj.Csym = obj.Parser.Csym;
            obj.Dsym = obj.Parser.Dsym;
        end
        
        
    end
end

