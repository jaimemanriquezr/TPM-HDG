classdef IdentityOperator < WeightOperator
    % WeightOperator with weight 1: applies no change to the value.
    methods
        function obj = IdentityOperator()
            obj = obj@WeightOperator(1.0);
        end
    end
end