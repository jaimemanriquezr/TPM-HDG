function new_obj = minus(obj, A)
arguments
    obj
    A
end
if ~isa(obj, "SimplicialGrid")
    [obj, A] = deal(A, obj);
end
new_obj = obj;
new_obj.coordinates = obj.coordinates - A;
end
