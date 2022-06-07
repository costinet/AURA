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
        
        circuitParser
        swseq
        
        filename % The filename that contains the topology for this covnerter
        
        stateLabels={}
        stateLabels_Opp={}
        Xi
        Bi
        outputLabels
        
        
        ts
        order % Need right now becuase it breaks code if not here
        Element_Properties % 1st column is char of all element names 2nd column is the value of that element
        Switch_Resistors  % List of chars that represent the switch names plus '_R'
        Switch_Resistor_Values % [SW_OFF; SW_ON; SW_ON] third one will eventually be diode resistance
        Switch_Sequence % Set of binary on or off values that has the number of colums of swithces and the numer of rows of time intervals
        Fwd_Voltage
        
    end
    
    % %%From DC Code
    % properties
    %     As
    %     Bs
    %     Cs
    %     Ds
    %     Is % Matrices denoting any dependent states
    %     K % Matrix of state component values (L's and C's)
    %
    %     swseq
    %
    %     circuitParser
    % end
    %
    % properties (Dependent)
    %     constraints
    %     labels
    %
    %         sourcefn
    %     end
    %
    %     properties (Hidden)
    %        % Constraints
    %         Cbnd
    %         Dbnd
    %         bndHyst     %[mean, std]
    %         switchRef   %[SwitchNo, State its in]
    %
    %         stateLabels
    %         outputLabels
    %         inputLabels
    %         switchLabels
    %     end
    %
    % %     properties (SetAccess = private)
    % %         plecsfn
    % %         plecsfdate
    % %     end
    %
    %     properties (Hidden, SetAccess = private)
    %         expAs
    %         eigAs
    %     end
    
    
    
    
    methods
        loadModelLT(obj, fn, swseq, force)
        % From external files
        function loadCircuit(obj, fn, swvec, force)
            
%             if nargin > 2
%                 obj.swseq = [obj.swseq; setdiff(swvec, obj.swseq, 'rows')];
%             end
            if nargin < 4
                force = 0;
            end
            
            if isempty(fn) && ~isempty(obj.circuitParser) && ~isempty(obj.sourcefn)
                fn = obj.sourcefn;
            end
            
            if isempty(obj.circuitParser)
                try %Test if the file is a valid plecs circuit
                    plecs('get', fn, 'StateSpaceOrder');
                    obj.circuitParser = PLECScircuitParser(obj);
                catch e
                end
                
                if isempty(obj.circuitParser)
                    parse = NetListParse();
                    obj.circuitParser = parse;
                    parse.initialize(fn,{},{});
                    parse.cutset_loop_num();
                    obj.loadModelLT(fn, swvec, force);
                end
            end
            loadModel(obj.circuitParser, fn, swvec, force);
            obj.circuitParser.getConstraintMatrices;
        end
        
        
        
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

