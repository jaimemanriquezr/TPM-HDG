function obj_plus_x = plus(obj, x)
%PLUS  Sum of an FE expression and another expression or datum.
%
%   A double or function handle operand is shifted in as a Datum. Two
%   expressions with the same role (or both trial-bearing) form a
%   LinearCombination with coefficients [1 1]; otherwise the sum goes
%   through AbstractFEFunction.apply.
    arguments
        obj
        x
    end
    if ~isa(obj, 'AbstractFEFunction')
        x_ = obj;
        obj = x;
        x = x_;
    end

    switch class(x)
        case {'double', 'function_handle'}
            obj_plus_x = FEFunction(@(varargin) ...
                obj.evaluate(varargin{:}) + Datum(x).evaluate(varargin{:}), ...
                obj.signature);
        case "LinearCombination"
            obj_plus_x = plus(x, obj);
        otherwise
            if isequal(extract(obj.signature, lettersPattern), extract(x.signature, lettersPattern))
                obj_plus_x = LinearCombination([obj x], [1 1]);
                obj_plus_x.signature = join(extract(obj.signature, lettersPattern), ":");
            elseif all(contains([obj.signature, x.signature], "Trial"))
                % trial-bearing expressions with different wrappings still
                % form an inspectable sum (multi-field path_gradient, ...);
                % form signatures never contain "Trial", so this cannot
                % capture LinearForm/BilinearForm sums
                obj_plus_x = LinearCombination([obj x], [1 1]);
                obj_plus_x.signature = AbstractFEFunction.match_signature(obj, x);
            else
                obj_plus_x = AbstractFEFunction.apply(@plus, obj, x);
            end
    end
end