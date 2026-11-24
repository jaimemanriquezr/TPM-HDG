function obj_t = ctranspose(obj)
%CTRANSPOSE  Page-wise conjugate transpose; the id becomes [test trial].
    fun = @pagectranspose;
    array = fun(double(obj));
    if length(size(array)) > 5
        array = permute(array, [1 2 3 4 5 7 6]);
    end
    obj_t = IndexedArray(array, fliplr(obj.id));
end