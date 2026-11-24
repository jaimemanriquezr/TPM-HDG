function a_minus_b = minus(a, b)
%MINUS  Difference of an IndexedArray and another operand.
%
%   Arrays of equal id subtract to an IndexedArray; arrays of different
%   ids give an IndexedMatrix; b = 0 returns a; an IndexedMatrix b is
%   handled by its own minus; any other b gives a plain array.
    arguments
        a IndexedArray
        b
    end
    fun = @minus;
    if isa(b, 'IndexedMatrix')
        a_minus_b = -minus(b, a);
        return
    end

    if isequal(b, 0.0)
        a_minus_b = a;
    elseif isa(b, 'IndexedArray')
        if isequal(a.id, b.id)
            a_minus_b = IndexedArray(fun(double(a), double(b)), a.id);
        else
            a_minus_b = fun(IndexedMatrix(a), IndexedMatrix(b));
        end
    else
        a_minus_b = fun(double(a), double(b));
    end
end