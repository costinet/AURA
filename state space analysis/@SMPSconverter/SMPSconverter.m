classdef SMPSconverter < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        topology
        ts
        u
        
        swind %indices into the A/B/C/D etc. matrices stored by the topology that 
        %form the switching pattern. e.g. swseq = [1 3] will form a period using 
        %As(:,:,1) and As(:,:,3)
        
        tcomp % compensation indexes. When ts(i) decreases by dt, ts(tcomp(i)) should increase by dt
        tsmax % max permissible values for ts
    end
    
    properties (Dependent = true)
        As
        Bs
        Cs
        Ds
        Is
        swseq = 1;
    end
    
    properties (SetAccess = private)
        expAts
        cachets
        cacheAs
    end
    
    methods
        function [tnew] = adjustTiming(obj, ti, dt)
            tnew = obj.ts;
            
            if( dt + tnew(ti) < 0)
                dt = -tnew(ti);
            end
            if( dt + tnew(ti) > obj.tsmax(ti))
                dt = obj.tsmax(ti) - tnew(ti);
            end
            if( tnew(obj.tcomp(ti)) - dt > obj.tsmax(obj.tcomp(ti)))
                dt = -obj.tsmax(obj.tcomp(ti)) + tnew(obj.tcomp(ti));
            end
           
            tnew(ti) = tnew(ti) + dt;
            tnew(obj.tcomp(ti)) = tnew(obj.tcomp(ti)) - dt;
            
            assert(~(sum(tnew < 0) && sum(tnew>obj.tsmax)), 'Timing adjustment failed');
            
            obj.ts = tnew;
            obj.refreshCache();
        end
        
        function refreshCache(obj)
            if(size(obj.swseq) == size(obj.As,3))
                for i = 1:size(obj.As,3)
                    obj.expAts(:,:,i) = expm(obj.As(:,:,obj.swseq(i))*obj.ts(i));
                end 
            end
        end
        
        %% Getters
        function res = get.As(obj)
            res = obj.topology.As(:,:,obj.swseq);
        end
        
        function res = get.Bs(obj)
            res = obj.topology.Bs(:,:,obj.swseq);
        end
        
        function res = get.Cs(obj)
            res = obj.topology.Cs(:,:,obj.swseq);
        end
        
        function res = get.Ds(obj)
            res = obj.topology.Ds(:,:,obj.swseq);
        end
        
        function res = get.Is(obj)
            res = obj.topology.Is(:,:,obj.swseq);
        end
        
        function res = get.swseq(obj)
            if isempty(obj.swind)
                res = 1:size(obj.topology.As,3);
            else
                res = obj.swind;
            end
        end
        
        function expAts = get.expAts(obj)
            expAts = obj.expAts;
        end
        
        %% Setters
        function set.swseq(obj,newSeq)
            obj.swind = newSeq;
            obj.refreshCache();
        end
        
        function set.ts(obj, newts)
            obj.ts = newts;
            obj.refreshCache();
        end
    end
    
end

