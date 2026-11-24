classdef TracedFunction < AbstractFEFunction
    % trace(u) as a structured node (not a closure): keeps the operand
    % walkable for expression walkers while evaluating identically.
    properties
        u
    end
    methods
        function obj = TracedFunction(u)
            obj.u = u;
            obj.signature = u.signature;
        end
        function p = evaluate(obj, varargin)
            p = TraceOperator.apply_action(evaluate(obj.u, varargin{:}));
        end
        function p = evaluate_boundary(obj, varargin)
            p = TraceOperator.apply_action(evaluate_boundary(obj.u, varargin{:}));
        end
    end
end
