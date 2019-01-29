classdef FET < handle
    %FET contains all of the data of a Field Effect Transistor (FET)
    %   Still needs some work
    
%     %%%%%%   %      %  %%%%%%%    %%%%%%
%    %      %  %      %  %      %  %      %
%    %      %  %      %  %      %  %      %
%    %%%%%%%%  %      %  %%%%%%%   %%%%%%%%
%    %      %  %      %  %%        %      %
%    %      %  %      %  % %       %      %
%    %      %  %      %  %  %      %      %
%    %      %  %      %  %   %     %      %
%    %      %   %    %   %    %    %      %
%    %      %    %%%%    %     %   %      %
    
    properties

        R_on
        R_diode
        V_f
        R_off
        
    end

    methods

        function initialize(obj,Filename,Voltage,Current)
            % This function checks to see if all variables are of the
            % correct class
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
