classdef Results < handle
    %Results Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Xs
        ts
        ONorOFF
    end
    
    methods
        
        function [pass]=checkresults(obj,Xs,ts,ONorOFF)
            
             XsCheck = norm((obj.Xs-Xs)./obj.Xs);
             tsCheck = norm((obj.ts-ts)./obj.ts);
             ONorOFFCheck = sum(sum(obj.ONorOFF-ONorOFF));
            
             if XsCheck+tsCheck+ONorOFFCheck<0.1
                 pass = 1;
             else 
                 pass = 0;
             end
                 
        end
        
    end
    
end

