classdef component < smps.databases.component
     % Shim class for loading old AURAdb data
    
    methods (Static)
        function obj = loadobj(A)
            warning('off', 'MATLAB:load:classDoesNotMatch');
            obj = smps.components.component();
            obj.partNumber = A.partNumber;
            obj.merge(A);
            warning('on', 'MATLAB:load:classDoesNotMatch');
        end
    end

end