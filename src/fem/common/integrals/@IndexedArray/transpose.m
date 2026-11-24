function obj_t = transpose(obj)
%TRANSPOSE  Page-wise transpose; the id becomes [test trial].
    fun = @pagetranspose;
    array = fun(double(obj));
    if length(size(array)) > 5
        array = permute(array, [1 2 3 4 5 7 6]);
    end
    obj_t = IndexedArray(array, fliplr(obj.id));
end