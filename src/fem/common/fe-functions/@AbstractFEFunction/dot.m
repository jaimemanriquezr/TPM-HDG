function obj_dot_x = dot(obj, x)
%DOT  Inner product dot(u, x) of an FE expression with another operand.
%
%   x a WeightFunction gives a Product; x a numeric array contracts the
%   value with it (DotOperator); any other x goes through
%   AbstractFEFunction.apply (a BilinearForm for a trial/test pair, a
%   LinearForm for a test against data). Either operand may be the FE
%   expression.
    arguments
        obj
        x
    end
    if ~isa(obj, "AbstractFEFunction")
        x_ = obj;
        obj = x;
        x = x_;
    end
    
    if isa(x, "WeightFunction")
        obj_dot_x = Product([obj, x], @AbstractFEFunction.dot_);
        obj_dot_x.signature = obj.signature;
    else
        switch class(x)
            case {"double"}
                obj_dot_x = FEFunction(@(varargin) ...
                DotOperator(x).evaluate(obj, varargin{:}), obj.signature);
            otherwise
                obj_dot_x = AbstractFEFunction.apply(@dot, obj, x);
        end
    end
end
