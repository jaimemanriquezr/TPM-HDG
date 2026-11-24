function w_u = times(w, u)
%TIMES  Elementwise product w.*u of a weight with a double or an FE expression.
%
%   A double u gives an array of scaled weights (one per entry of u);
%   anything else is multiplied as times(u, w).
    if ~isa(w, "WeightFunction")
        [w, u] = deal(u, w);
    end

    switch class(u)
        case "double"
            w_u = arrayfun(@(u) mtimes(u, w), u);
        otherwise
            w_u = times(u, w);
    end
end