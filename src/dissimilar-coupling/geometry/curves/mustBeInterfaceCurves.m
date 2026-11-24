function mustBeInterfaceCurves(interface)
%MUSTBEINTERFACECURVES  One AbstractCurve, or a pair {curve_1, curve_2} on one parameter.
%
%   Argument validator of the TPM couplers: a scalar AbstractCurve, or a
%   1 x 2 cell of scalar AbstractCurves with equal t_range and equal
%   periodicity (see interface_curves).
if isa(interface, "AbstractCurve") && isscalar(interface)
    return
end
pair = iscell(interface) && numel(interface) == 2 ...
    && all(cellfun(@(c) isa(c, "AbstractCurve") && isscalar(c), interface));
if ~pair
    error("mustBeInterfaceCurves:type", ...
        "The interface must be an AbstractCurve or a pair {curve_1, curve_2} of them.");
end
if ~isequal(interface{1}.t_range, interface{2}.t_range) ...
        || is_periodic(interface{1}) ~= is_periodic(interface{2})
    error("mustBeInterfaceCurves:parameter", ...
        "The two curves of a pair must share the parameter domain " + ...
        "(t_range and periodicity).");
end
end
