classdef topologyDB < smps.databases.topologyDB
    %Shim class for loading old AURAdb data

    methods (Static)
        function obj = loadobj(A)
            warning('off', 'MATLAB:load:classDoesNotMatch');
            obj = smps.databases.topologyDB(1);
            obj.addMult(A.components)
            warning('on', 'MATLAB:load:classDoesNotMatch');
        end
    end

end