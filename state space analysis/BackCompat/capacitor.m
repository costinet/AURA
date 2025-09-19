classdef (Hidden) capacitor < smps.components.capacitor
     % Shim class for loading old AURAdb data
    
    methods (Static)
        function c = loadobj(A)
            warning('off', 'MATLAB:load:classDoesNotMatch');
            c = smps.components.capacitor();
            c.partNumber = A.partNumber;
            c.merge(A);
            warning('on', 'MATLAB:load:classDoesNotMatch');
        end
    end

end