classdef FEOperator < matlab.mixin.Heterogeneous
    % Pointwise operator acting on the evaluated value of an FE function:
    % action(p), or action(p, x) when it depends on the position
    % (is_spatial_dependent). input_dimension and output_dimension are the
    % value sizes, -1 when unspecified.
    properties
        action function_handle
        input_dimension double = -1;
        output_dimension double = -1;
    end

    properties (Dependent)
        is_spatial_dependent
    end

    methods
        function obj = FEOperator(action)
            arguments
                action function_handle = @(p) p
            end
            obj.action = action;
        end
        
        function bool = get.is_spatial_dependent(obj)
            bool = (nargin(obj.action) > 1);
        end
        
        function F_p = apply_action(obj, p, varargin)
            F_p = obj.action(p);
        end

        F_p = evaluate(obj, u, x_local, mesh, options)
    end
end
