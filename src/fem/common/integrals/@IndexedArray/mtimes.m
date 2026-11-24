function x = mtimes(a, b)
%MTIMES  Product b*double(a) with the IndexedArray as the second factor.
%
%   The operands are swapped if the first is not an IndexedArray. For a
%   double b the result is an IndexedArray with the id of a; otherwise a
%   plain array.
    if ~isa(a, "IndexedArray")
        [a, b] = deal(b, a);
    end
    
    if isa(b, "double")
        x = IndexedArray(b * double(a), a.id);
    else
        x = b * double(a);
    end
end