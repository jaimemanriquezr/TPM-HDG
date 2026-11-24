function obj = subsasgn(obj, s, b)
%SUBSASGN  Indexed assignment on the array part; the id is kept.
    if isequal(s.type, '.')
        obj = builtin('subsasgn', obj, s, b);
    else
        obj = IndexedArray(builtin('subsasgn', double(obj), s, b), obj.id);
    end
end