function s = sum(obj)
    % contract an elementwise product: sum(S.*L) == dot(S, L)
    arguments
        obj BilinearForm
    end
    assert(isequal(func2str(obj.op), "times"), "AbstractFEFunction:sum", ...
        "sum of a bilinear form needs an elementwise product (u.*v)");
    s = BilinearForm(obj.factors(1), obj.factors(2), @dot);
end
