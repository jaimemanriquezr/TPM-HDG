function b = subsref(obj, s)
%SUBSREF  Indexing of the array part, returning an IndexedArray with the same id.
%
%   Dot references (the id property) are passed to the builtin.
    if isequal(s.type, '.')
        b = builtin('subsref', obj, s);
    else
        b = IndexedArray(builtin('subsref', double(obj), s), obj.id);
    end
end