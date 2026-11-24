classdef IndexedArray < double
    % A numeric array (a double subclass) tagged with the id [trial test]
    % of the pair of spaces it couples. Arithmetic keeps the id; arrays of
    % different ids combine into an IndexedMatrix.
    properties
        id (1, 2) double
    end
    
    methods 
        function obj = IndexedArray(array, id)
            arguments
                array double = []
                id (1, 2) double = [0 0]
            end
            obj = obj@double(array);
            obj.id = id;
        end

        function obj = assemble(obj)
        end
    end
end