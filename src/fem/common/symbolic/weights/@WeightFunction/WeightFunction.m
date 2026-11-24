classdef WeightFunction < AbstractFEFunction
    % Base class of pointwise weights (data that multiplies FE functions).
    % The weight property holds the weight; a concrete subclass implements
    % evaluate, and evaluate_boundary defaults to it.
    properties (SetAccess=protected)
        weight
    end
    
    methods
        function w = evaluate(obj, x_local, mesh, options)
            arguments
                obj WeightFunction
                x_local double
                mesh AbstractGrid
                options.coordinate_system = "";
                options.output = "";
            end
            error("WeightFunction:abstract", ...
                "WeightFunction is a base class; evaluate a concrete subclass (%s).", class(obj));
        end
        
        function w = evaluate_boundary(obj, varargin)
            w = evaluate(obj, varargin{:});
        end
    end
end