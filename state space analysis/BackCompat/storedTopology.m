classdef storedTopology < smps.components.storedTopology
    % Shim class for loading old AURAdb data
    
    methods (Static)
        function t = loadobj(A)
            warning('off', 'MATLAB:load:classDoesNotMatch');
            t = smps.components.storedTopology();
            t.loadshim(A);
            warning('on', 'MATLAB:load:classDoesNotMatch');
        end
    end
end