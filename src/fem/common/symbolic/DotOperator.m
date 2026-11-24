classdef DotOperator < WeightOperator
    % WeightOperator that contracts the value of its argument with the
    % weights (a numeric array or a weight function evaluated at the
    % physical points).
    methods
        function obj = DotOperator(weights)
            arguments
                weights = 0
            end
            obj = obj@WeightOperator(weights);
        end
        
        function W_p = apply_action(obj, p, x)
            arguments
                obj WeightOperator
                p double
                x double = []
            end
            v_size = size(p, 3:4);
            if isa(obj.weights, 'double')
                w = reshape(obj.weights, v_size);
                w_p_ = sum(w .* permute(p, [3 4 1 2 5 6 7]), 1:2);
                W_p = permute(w_p_, [3 4 1 2 5 6 7]);
            else
                x_weights = obj.weights.evaluate(x, coordinate_system = "physical");
                w = reshape(x_weights, [size(x_weights, 1:2) v_size size(x_weights, 5)]);
                w_p_ = sum(w .* permute(p, [3 4 1 2 5 6 7]), 1:2);
                W_p = permute(w_p_, [3 4 1 2 5 6 7]);
            end
        end
    end
end