classdef (Hidden) capacitorDB < smps.databases.capacitorDB
    %Shim class for loading old AURAdb data

    methods (Static)
        function obj = loadobj(A)
            warning('off', 'MATLAB:load:classDoesNotMatch');
            obj = smps.databases.capacitorDB();
            obj.addMult(A.components)
            warning('on', 'MATLAB:load:classDoesNotMatch');
        end
    end

end