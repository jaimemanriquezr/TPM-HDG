function obj_minus_x = minus(obj, x)
%MINUS  Difference of an FE expression and another expression or datum.
%
%   A double or function handle operand is shifted in as a Datum. Two
%   expressions with the same role (or both trial-bearing) form a
%   LinearCombination with coefficients [1 -1]; otherwise the difference
%   goes through AbstractFEFunction.apply.
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
            obj_minus_x = FEFunction(@(varargin) ...
                obj.evaluate(varargin{:}) - Datum(x).evaluate(varargin{:}), ...
                obj.signature);
        case "LinearCombination"
            % obj - x with x a LinearCombination (minus(x, obj) would
            % compute x - obj)
            obj_minus_x = plus(uminus(x), obj);
        otherwise
            if isequal(extract(obj.signature, lettersPattern), extract(x.signature, lettersPattern))
                obj_minus_x = LinearCombination([obj x], [1 -1]);
                obj_minus_x.signature = join(extract(obj.signature, lettersPattern), ":");
            elseif all(contains([obj.signature, x.signature], "Trial"))
                % see plus: inspectable sum for trial-bearing expressions
                obj_minus_x = LinearCombination([obj x], [1 -1]);
                obj_minus_x.signature = AbstractFEFunction.match_signature(obj, x);
            else
                obj_minus_x = AbstractFEFunction.apply(@minus, obj, x);
            end
    end
end