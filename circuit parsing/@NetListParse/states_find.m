function [ON_States_ALL] = states_find(obj,state,switches)
% STATES finds the different switch positions of the netlist file converter
% and sets the net list for nodeloop by changing open switches to caps and
% on switches to caps in parallel with an on resistance.


%     _   _   _  ____    _
%    / \ | | | |/ _  |  / \
%   / _ \| | | | (_| | / _ \
%  / ___ | |_| |> _  |/ ___ \
% /_/   \_\___//_/ |_/_/   \_\


% Does the current implementation have a forward voltage state and a
% non forward voltage state for diodes????

% Index values for components (will pass as variable... eventually):
numV = 1;
numBV = 2;
numMV = 3;
numC = 4;
numR = 5;
numL = 6;
numMI = 7;
numBI = 8;
numI = 9;
numD = 10;
numM = 11;

NL = obj.NL;
NLnets = obj.NLnets;
ON_States = cell(size(switches,2),1);
ON_index = 1;



for i = 1:1:size(state,1)
    ON_States = cell(size(switches,2),1);
    ON_index = 1;
    
    for j = 1:1:length(switches)
        
        % OPEN switch is replaced with cap
        if j~=length(switches)
            
            % if a FET then need to check that state and next state is
            % zero AND check that the current and next states correspond
            % to the same switch AND that is is a FET
            
            
            if (state(i,j) == 0 && state(i,j+1) == 0 && switches(j+1) == switches(j) && (NL(switches(j),1)==numM))

            % Add resistor in parallel to Cap:
            % Add new voltage measurement element to NLnets:
            NLnets_State = [strcat(NLnets(switches(j),1),'_R_OFF')];
            
            ON_States(ON_index,1) = NLnets_State;
            ON_index = ON_index+1;
            end
            
            % CLOSED switch is replaced with cap and resistor
            if (state(i,j) == 1 && switches(j+1) == switches(j)) && (NL(switches(j),1)==numM)

            
            % Add resistor in parallel to Cap:
            % Add new voltage measurement element to NLnets:
            NLnets_State = [strcat(NLnets(switches(j),1),'_R_ON')];
            
            ON_States(ON_index,1) = NLnets_State;
            ON_index = ON_index+1;
                
                % BODY DIODE switch is replaced with resistor and forward voltage
                % source have different if statements for diodes and FETs due to
                % FETs diode being opposite in direction given node and sign
                % convention for LTSpice
            end
            
        end
        %{
        % if Diode then need to check that it is a diode and that the
        % state is set to zero
        if state(i,j) == 0 && (NL(switches(j),1)==numD)

            % Add resistor in parallel to Cap:
            % Add new voltage measurement element to NLnets:
            NLnets_State = [strcat(NLnets(switches(j),1),'_R_OFF')];
            
            ON_States(ON_index,1) = NLnets_State;
            ON_index = ON_index+1;
            
        end
        
        
        % For power diodes
        
        % Check if state is set to 1 and that it is a diode
        if state(i,j) == 1 && (NL(switches(j),1)==numD)

            % Add resistor in parallel to Cap:
            % Add new voltage measurement element to NLnets:
            NLnets_State = [strcat(NLnets(switches(j),1),'_R_D')];
            
            ON_States(ON_index,1) = NLnets_State;
            ON_index = ON_index+1;
            
            
        end
        
        % For FETS body diode ON
        if j>1 % Without this have an index error because need to check the previous index
            % Have to check if the current state is set to 1 AND if the
            % previous index represents the same state AND if the index
            % corresponds to a FET
            if state(i,j) == 1 && (switches(j-1) == switches(j)) && (NL(switches(j),1)==numM)

                % Add resistor in parallel to Cap:
                % Add new voltage measurement element to NLnets:
                NLnets_State = [strcat(NLnets(switches(j),1),'_R_D')];
                
                ON_States(ON_index,1) = NLnets_State;
                ON_index = ON_index+1;
                
            end
        end
        
        %}
        if state(i,j) ~= 1 && state(i,j) ~= 0
            error('Invalid State Detected')
        end
        
    end
    
    ON_States_ALL(:,i)=ON_States;
end



end % That's all Folks
