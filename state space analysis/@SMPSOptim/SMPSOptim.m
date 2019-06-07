classdef SMPSOptim < handle
    %OPTIM contains all of the data for the optimization of the
    %converter
    %   Just getting started
    



%     _   _   _  ____    _    
%    / \ | | | |/ _  |  / \   
%   / _ \| | | | (_| | / _ \  
%  / ___ | |_| |> _  |/ ___ \ 
% /_/   \_\___//_/ |_/_/   \_\




properties
    
    ts_no_dead
    ts
    
    Simulator
    
    dead_time_intervals
    dead_time_states
    dead_time_goals
    
    Vo_index
    Vo_ideal_value
    Perturb1_index
    Perturb2_index
    
end

    methods
        [dXs] = fs_adjust(obj)
        function initialize(obj)
            % This function creates variables for all of the FET
            % properties
            
            
            
            
            
            if ~exist(Filename,'file')
                error('File %s not found',Filename)
            end

            if ischar(Filename)
                obj.filename=Filename;
            else
                error('Filename must be of type char \nCurrent class of filename: %s',class(Filename))
            end
            if nargin == 2
                obj.Meas_Voltage = {};
                obj.Meas_Current = {};
                return
            end
            if nargin == 3
                warning('Only voltage measurements detected')
                Current = {};
            end
            if iscell(Voltage)
                if size(Voltage,1)==1||size(Voltage,2)==1||isempty(Voltage)
                    obj.Meas_Voltage = Voltage;
                else
                    error('Voltage is not the correct dimension')
                end
            else
                error('Voltage must be of type cell \nCurrent class of Voltage: %s',class(Voltage))
            end
            if iscell(Current)
                if size(Current,1)==1||size(Current,2)==1||isempty(Current)
                    obj.Meas_Current = Current;
                else
                    error('Current is not the correct dimension')
                end
            else
                error('Current must be of type cell \nCurrent class of Current: %s',class(Current))

            end
        end

        function [StateNames]=getStateNames(obj)
            StateNames = obj.StateNames;
        end


    end
    % That's all Folks
end
