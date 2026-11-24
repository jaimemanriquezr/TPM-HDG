function [fluid_mesh, solid_mesh] = fsi_dsl_meshes(interface, family, h, options)
% The (fluid, solid) mesh pair of the steady FSI drivers; the fluid lies
% on the upper/outer side of the interface, the solid on the lower/inner.
%
%   "fitted"  both sides meshed to the linear interpolant of Sigma
%             (mesh size h); gap > 0 retreats both polygons from it
%             (vertically for a graph, along the normal for a closed
%             curve).
%   "active"  active meshes of two crossed background grids of the unit
%             square (each a subset of its grid): the elements Sigma cuts
%             are discarded, so a gap band holds Sigma. The solid grid has
%             size h, the fluid grid h/2.
arguments
    interface AbstractCurve
    family (1,1) string {mustBeMember(family, ["fitted", "active"])}
    h (1,1) double {mustBePositive}
    options.gap (1,1) double {mustBeNonnegative} = 0
end
switch family
    case "fitted"
        [fluid_mesh, solid_mesh] = interface_meshes(interface, "fitted", h, ...
            gap = options.gap);
    case "active"
        assert(options.gap == 0, "fsi_dsl_meshes:gap", ...
            "Active meshes take no gap.");
        N = round(1/h);
        bg_fluid = StructuredSimplicialGrid(2*N, 2*N, "crossed");
        bg_solid = StructuredSimplicialGrid(N, N, "crossed");
        [fluid_mesh, ~] = bg_fluid.split(interface, cut_elements = "discard");
        [~, solid_mesh] = bg_solid.split(interface, cut_elements = "discard");
end
end
