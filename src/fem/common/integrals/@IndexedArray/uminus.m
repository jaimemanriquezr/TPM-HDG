function minus_a = uminus(a)
%UMINUS  Negation; the id is kept.
    arguments
        a IndexedArray
    end
    fun = @uminus;
    minus_a = IndexedArray(fun(double(a)), a.id);
end