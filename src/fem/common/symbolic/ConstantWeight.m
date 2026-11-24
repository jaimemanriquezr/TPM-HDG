classdef ConstantWeight < WeightFunction
    % A constant scalar/matrix weight as a structured node: scalar
    % coefficients and constant matrix actions (c*u, kappa\u) build
    % Product([ConstantWeight(A), u], ...) instead of an opaque closure,
    % so expression walkers (the transmission forms) can still see the
    % underlying primitive. Value shape rides in dims 3-4.
    methods
        function obj = ConstantWeight(A)
            obj.weight = A;
            obj.signature = "ConstantWeight";
        end
        function w = evaluate(obj, varargin)
            A = obj.weight;
            w = reshape(A, [1, 1, size(A, 1), size(A, 2)]);
        end
        function w = evaluate_boundary(obj, varargin)
            w = evaluate(obj);
        end
    end
end
