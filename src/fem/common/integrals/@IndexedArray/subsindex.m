function ind = subsindex(A)
%SUBSINDEX  Use of an IndexedArray as an index: its numeric values.
    ind = builtin('subsindex', double(A));
end