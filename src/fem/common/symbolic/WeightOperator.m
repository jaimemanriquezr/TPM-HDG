classdef WeightOperator < FEOperator
    % Pointwise weighting operator: (W u)(x) = w(x) * u(x).
    %
    % Weight forms:
    %   double        constant scalar or matrix, e.g. WeightOperator(inv(kappa))
    %   spatial       anything Datum accepts: a scalar function handle
    %                 @(x, y), a sym, or a cell matrix of scalar handles for
    %                 a matrix-valued weight, e.g.
    %                 WeightOperator({@(x,y) 1+x.^2, @(x,y) 0*x; ...})
    %
    % A matrix weight applied to a scalar-valued operand is promoted to
    % elementwise scaling (p*[0 1; -1 0], trace(sigma)*eye(2), ...);
    % otherwise weights compose by matrix multiplication on the value.
    %
    % Quadrature rules are chosen from the polynomial orders of the bases
    % alone; a non-polynomial spatial weight is integrated inexactly at the
    % default orders. Raise the rule degree at the call site when needed.
    properties
        weights
    end

    methods
        function obj = WeightOperator(weights)
            arguments
                weights = 0
            end
            % The closures are shadowed by apply_action, but their arity
            % matters: FEOperator.is_spatial_dependent tests nargin(action)
            % > 1, which routes physical points into apply_action.
            if isa(weights, 'double')
                id_action = @(p) weights * p;
            else
                weights = Datum(weights);
                id_action = @(p, x) weights.evaluate(x) * p;
            end
            obj = obj@FEOperator(id_action);
            obj.weights = weights;
        end

        function W_p = apply_action(obj, p, x)
            arguments
                obj WeightOperator
                p double
                x double = []
            end
            p_ = permute(p, [3 4 1 2 5 6 7]);   % (o1 x o2 x N x d x ...)
            if isa(obj.weights, 'double')
                w = obj.weights;
            else
                w = obj.weights.evaluate(x, coordinate_system = "physical");
                w = permute(w, [3 4 1 2 5 6 7]); % (s1 x s2 x N x 1 x ...)
            end
            W_p = permute(WeightOperator.compose(w, p_), [3 4 1 2 5 6 7]);
        end
    end

    methods (Static)
        function w_p = compose(w, p_)
            % w and p_ carry values in dims 1:2 and points/dofs in the
            % page dims; scalar values scale, matching shapes multiply
            if size(p_, 1) == 1 && size(p_, 2) == 1 || isscalar(w(:, :, 1))
                w_p = w .* p_;
            else
                w_p = pagemtimes(w, p_);
            end
        end
    end
end
