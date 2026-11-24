function [curve_1, curve_2] = interface_curves(interface)
%INTERFACE_CURVES  The per-side curves of a TPM interface.
%
%   [curve_1, curve_2] = interface_curves(interface)
%
%   interface is one AbstractCurve, shared by both sides (curve_2 is then
%   []), or a pair {curve_1, curve_2}: side 1 maps its mesh onto curve_1
%   and side 2 onto curve_2, and the two are paired through their common
%   parameter t (get_partitions, option interface_2). The pair must share
%   t_range and periodicity. For the fsi coupling the pair is
%   {curve_fluid, curve_solid}: the fluid sees the current interface
%   x(t) and the solid its reference X(t) (see InterfaceProblem.fsi).
%
%See also: mustBeInterfaceCurves, AbstractCurve.get_partitions
mustBeInterfaceCurves(interface);
if isa(interface, "AbstractCurve")
    curve_1 = interface;
    curve_2 = [];
else
    [curve_1, curve_2] = interface{:};
end
end
