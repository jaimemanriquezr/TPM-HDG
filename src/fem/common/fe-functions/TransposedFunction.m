classdef TransposedFunction < AbstractFEFunction
    % u' or u.' as a structured node (not a closure): keeps the operand's
    % evaluate_boundary, so mesh-weighted products (tau*u, sigma*n)
    % transpose correctly on element boundaries.
    properties
        u
        op function_handle = @AbstractFEFunction.ctranspose_
    end
    methods
        function obj = TransposedFunction(u, op)
            obj.u = u;
            obj.op = op;
            obj.signature = u.signature;
        end
        function p = evaluate(obj, varargin)
            p = obj.op(evaluate(obj.u, varargin{:}));
        end
        function p = evaluate_boundary(obj, varargin)
            p = obj.op(evaluate_boundary(obj.u, varargin{:}));
        end
    end
end
