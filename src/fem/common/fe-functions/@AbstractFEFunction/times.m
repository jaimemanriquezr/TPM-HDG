function obj_times_x = times(obj, x)
    % elementwise product u.*v; contract with sum (sum(S.*L) == dot(S, L))
    arguments
        obj
        x
    end
    if ~isa(obj, "AbstractFEFunction")
        x_ = obj;
        obj = x;
        x = x_;
    end

    if isa(x, "double")
        if isscalar(x)
            obj_times_x = obj * x;
        else
            obj_times_x = Product([ConstantWeight(x), obj], @times);
            obj_times_x.signature = obj.signature;
        end
    else
        obj_times_x = AbstractFEFunction.apply(@times, obj, x);
    end
end
