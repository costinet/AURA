classdef transistor < smps.components.transistor
    % Shim class for loading old AURAdb data
    
    methods (Static)
        function t = loadobj(A)
            warning('off', 'MATLAB:load:classDoesNotMatch');
            t = smps.components.transistor();
            t.partNumber = A.partNumber;
            t.merge(A);
            warning('on', 'MATLAB:load:classDoesNotMatch');
        end
    end
end