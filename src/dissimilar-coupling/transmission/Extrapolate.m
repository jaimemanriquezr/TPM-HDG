function obj = Extrapolate(x)
% Extrapolate(u): evaluate the local trial expression at the interface
% points x + m(x) of its side's induced patch (the TPM extrapolation
% E(u)). Returns a ExtendedFunction; a plain constructor function so the
% single ExtendedFunction classdef carries the (non-inherited)
% InferiorClasses dominance list. See ExtendedFunction.
obj = ExtendedFunction(x);
for i = 1:numel(obj.terms)
    assert(obj.terms(i).kind == "onface", "Extrapolate:arg", ...
        "Extrapolate wraps local (volume) trials only");
    obj.terms(i).kind = "extrapolate";
end
end
