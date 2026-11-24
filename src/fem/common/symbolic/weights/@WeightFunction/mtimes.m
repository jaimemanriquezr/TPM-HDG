function w_u = mtimes(obj, u)
%MTIMES  Product of a weight with a datum, a combination or an FE function.
%
%   A double scales the element weights; a LinearCombination has its
%   coefficients scaled; an FEField gives a Product; anything else is
%   multiplied as mtimes(u, w). A weight of 1 returns u, and a form on the
%   right is contracted. The operands are swapped if the first is not a
%   WeightFunction.
    if isa(u, "BilinearForm")
        % a weight on the left of a form contracts its value (L*(a*b))
        w_u = contract(u, obj, "left");
        return
    end
    if ~isa(obj, "WeightFunction")
        [obj, u] = deal(u, obj);
    end
    
    if isequal(obj.weight, 1)
        w_u = u;
        return
    end

    switch class(u)
        case "double"
            w_u = obj;
            for i = 1:length(w_u)
                w_u(i).element_weights = u * obj(i).element_weights;
            end
        case "LinearCombination"
            w_u = LinearCombination(u.summands, obj .* u.coefficients);
            w_u.signature = u.signature;
        otherwise
            if contains(u.signature, "FEField")
                w_u = Product([obj u], @times);
                w_u.signature = u.signature;
            else
                w_u = mtimes(u, obj);
            end
    end
end