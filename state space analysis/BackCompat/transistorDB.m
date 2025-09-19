classdef transistorDB < smps.databases.transistorDB
    %Shim class for loading old AURAdb data

    methods (Static)
        function obj = loadobj(A)
            warning('off', 'MATLAB:load:classDoesNotMatch');
            obj = smps.databases.transistorDB();
            obj.addMult(A.components)
            warning('on', 'MATLAB:load:classDoesNotMatch');
        end
    end

end