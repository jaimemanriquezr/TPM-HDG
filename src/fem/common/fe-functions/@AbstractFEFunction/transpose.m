function obj_transpose = transpose(obj)
%TRANSPOSE  Transpose of the value dimensions, as a TransposedFunction.
    arguments
        obj AbstractFEFunction
    end
    obj_transpose = TransposedFunction(obj, @AbstractFEFunction.transpose_);
end
