function obj_ctranspose = ctranspose(obj)
%CTRANSPOSE  Conjugate transpose of the value dimensions, as a TransposedFunction.
    arguments
        obj AbstractFEFunction
    end
    obj_ctranspose = TransposedFunction(obj, @AbstractFEFunction.ctranspose_);
end
