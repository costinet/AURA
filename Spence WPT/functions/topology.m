classdef topology < handle
    %topology class to hold A B C and D

    properties
        As
        Bs
        Cs
        Ds
        Is
        K
        
        xNames
        yNames
        uNames
        swNames
        
        swseq

    end
    
    methods
        function setSS(obj, As, Bs, Cs, Ds, K)
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
            
            if(nargin>5)
                obj.K = K;
            else
                obj.K = eye(length(obj.As));
            end
        end
        

        
        function loadPLECsModel(obj, fn, swseq)
            ssOrder = plecs('get', fn, 'StateSpaceOrder');

            obj.swNames = cellfun(@(x) x(strfind(x,'FET')+3:end), ssOrder.Switches, 'un', 0)';
            obj.xNames = cellfun(@(x) x(strfind(x,'/')+1:end), ssOrder.States, 'un', 0);
            obj.yNames = cellfun(@(x) x(strfind(x,'/')+1:end), ssOrder.Outputs, 'un', 0);
            obj.uNames = cellfun(@(x) x(strfind(x,'/')+1:end), ssOrder.Inputs, 'un', 0);
            
            if nargin > 2
                obj.swseq = swseq;
                for i = 1:size(swseq,1)
                    plecs('set', fn, 'SwitchVector', swseq(i,:));
                    names = plecs('get', fn, 'Topology');
                    obj.As(:,:,i) = names.A;
                    obj.Bs(:,:,i) = names.B;
                    obj.Cs(:,:,i) = names.C;
                    obj.Ds(:,:,i) = names.D;
                    obj.Is(:,:,i) = names.I;
                end
            end

        end

              
    end
    
end

