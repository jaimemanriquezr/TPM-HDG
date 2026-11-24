function obj = TangentialProjection(x)
% TangentialProjection(u): apply the tangential projector (t x t) at the
% patch quadrature points, with t the rotated outward interface normal of
% each term's own side (the Beavers-Joseph-Saffman pairing). Accepts trace
% trials and Extrapolate/PathIntegral terms.
obj = ExtendedFunction.coerce(x);
for i = 1:numel(obj.terms)
    assert(obj.terms(i).normal == "", "TangentialProjection:arg", ...
        "term already carries a normal/tangent contraction");
    obj.terms(i).normal = "tangent";
end
end
