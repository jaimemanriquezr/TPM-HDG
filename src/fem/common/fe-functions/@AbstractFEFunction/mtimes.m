function obj_times_x = mtimes(obj, x)
%MTIMES  Product of an FE expression with a datum, weight, measure or form.
%
%   A BilinearForm on the right is contracted with the left operand; a
%   Measure integrates; a double becomes a ConstantWeight product (scalars
%   broadcast, matrices act on the value dimensions); a function handle is
%   a spatial datum; a WeightFunction forms a Product; any other operand
%   goes through AbstractFEFunction.apply (forms for trial/test pairs; a
%   scalar dot(a, b) of data times a trial u, in either order, is the
%   Product u*a.'*b, evaluated factor by factor on the element boundary).
    arguments
        obj
        x
    end
    if isa(x, "BilinearForm")
        % data on the left of a form contracts its value (L*(a*b))
        obj_times_x = contract(x, obj, "left");
        return
    end
    if ~isa(obj, "AbstractFEFunction")
        x_ = obj;
        obj = x;
        x = x_;
    end

    if isa(x, "Measure")
        obj_times_x = mtimes(x, obj);
    elseif isa(x, "double")
        if isequal(x, 1)
            obj_times_x = obj;
        else
            % structured node, not a closure: keeps the operand walkable
            % (transmission forms recover the primitive and its mesh).
            % Scalars broadcast (@times); matrices act on the left
            % (A*u contracting the value dims).
            if isscalar(x) || prod(value_size(obj)) == 1
                % scalars broadcast; a matrix weight on a scalar-valued
                % operand promotes elementwise (gamma*J, trace(sg)*eye(2))
                op = @times;
            else
                op = @AbstractFEFunction.mtimes_;
            end
            obj_times_x = Product([ConstantWeight(x), obj], op);
            obj_times_x.signature = obj.signature;
        end
    elseif isa(x, "function_handle")
        % a spatial datum, whatever the operand: a test operand gives a
        % linear form, any other an FE expression
        obj_times_x = AbstractFEFunction.apply(@mtimes, obj, Datum(x));
    elseif isa(x, "WeightFunction")
        % isa, not class(x): mesh functions (ExteriorNormal, ...) are
        % WeightFunction subclasses and must form a Product (q'*n)
        obj_times_x = Product([obj, x], @AbstractFEFunction.mtimes_);
        obj_times_x.signature = obj.signature;
    else
        obj_times_x = AbstractFEFunction.apply(@mtimes, obj, x);
    end
end