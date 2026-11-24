function a_plus_b = plus(a, b)
%PLUS  Sum of an IndexedArray and another operand.
%
%   Arrays of equal id add to an IndexedArray; arrays of different ids
%   give an IndexedMatrix; b = 0 returns a; an IndexedMatrix b is handled
%   by its own plus; any other b gives a plain array.
    arguments
        a IndexedArray
        b 
    end
    fun = @plus;
    if isa(b, 'IndexedMatrix')
        a_plus_b = fun(b, a);
        return
    end
    
    if isequal(b, 0.0)
        a_plus_b = a;
    elseif isa(b, 'IndexedArray')
        if isequal(a.id, b.id)
            a_plus_b = IndexedArray(fun(double(a), double(b)), a.id);
        else
            a_plus_b = fun(IndexedMatrix(a), IndexedMatrix(b));
        end
    else
        a_plus_b = fun(double(a), double(b));
    end
end
