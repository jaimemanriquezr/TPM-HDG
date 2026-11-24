function plus_a = uplus(a)
%UPLUS  Unary plus; the id is kept.
    arguments
        a IndexedArray
    end
    fun = @plus;
    plus_a = IndexedArray(fun(double(a)), a.id);
end