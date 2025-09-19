classdef componentTableData < smps.components.componentTableData
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    methods (Static)
        function obj = loadobj(A)
            obj = smps.components.componentTableData(A);
        end
    end

end