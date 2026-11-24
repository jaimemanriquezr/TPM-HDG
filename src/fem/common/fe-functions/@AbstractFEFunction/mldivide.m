function obj_div_x = mldivide(obj, x)
%MLDIVIDE  Left division of FE expressions.
%
%   A double or function handle divisor is applied as a multiplication by
%   inv(x); otherwise the page-wise mldivide goes through
%   AbstractFEFunction.apply. The operands are swapped when the first is
%   not an FE expression.
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
            obj_div_x = mtimes(obj, inv(x));
        otherwise
            obj_div_x = AbstractFEFunction.apply(@mldivide, obj, x);
    end
end
