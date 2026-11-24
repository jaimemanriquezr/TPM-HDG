function [fluid_mesh, solid_mesh, info] = interface_meshes(curve, family, h, options)
%INTERFACE_MESHES  The two meshes of a TPM-coupled problem, by family.
%
%   [fluid_mesh, solid_mesh] = interface_meshes(curve, family, h, ...)
%
%   One mesh construction for every interface geometry. The families:
%
%     fitted          both subdomains meshed to the same linear interpolant
%                     of Sigma, gap `gap` (0 = they share its vertices)
%     interior        a background grid split by Sigma with the cut
%                     elements discarded: each mesh lies strictly inside
%                     its own subdomain, a gap band holds Sigma
%     protrude_fluid  the same split, cut elements to the fluid side
%     protrude_solid  the same split, cut elements to the solid side
%
%   `fluid_side` says which side of Sigma is the fluid, so the returned
%   pair is always (fluid, solid) whatever the geometry: "inner"/"outer"
%   for a closed curve, "lower"/"upper" for a graph. The default, outer,
%   is load_problems' default assignment (swap_domains = false). The unfitted
%   families take a crossed StructuredSimplicialGrid on `box`, or, when
%   `size_function` is set, an unstructured mesh2d triangulation of `box`
%   graded by it (hfun(PP), PP N-by-2, as refine2 accepts); the carve is
%   the same either way, so the meshes stay unfitted to Sigma.
%
%   max_normal_deviation is left to split's own policy unless set;
%   forcing 45 on a protrude family erodes both meshes and tears them
%   apart. `info` carries the invariants that catch such a tear:
%   n_shared (interface faces common to both meshes) and area_gap
%   (|box| - |fluid| - |solid|), zero for fitted with gap 0 and for the
%   protrude families, ~2*gap*perimeter for a gapped fitted pair, and
%   O(h) for interior. mesh2d may refine a thin polygon's boundary on
%   one side only, so a fitted pair's n_shared can fall short of its
%   face count.
%
%   The naive "Sigma = the mesh's own facets" comparison is an interface
%   mapping, not a mesh family: see interface_polygon.
%
%See also: fitted_split, SimplicialGrid.split, interface_polygon
arguments
    curve  AbstractCurve
    family (1,1) string {mustBeMember(family, ...
        ["fitted", "interior", "protrude_fluid", "protrude_solid"])}
    h      (1,1) double {mustBePositive}
    options.fluid_side (1,1) string ...
        {mustBeMember(options.fluid_side, ["inner", "outer", "lower", "upper"])} = "outer"
    options.box (1,4) double = [0 1 0 1]
    options.gap (1,1) double = 0          % fitted only: retreat from Sigma
    % fitted only: the solid mesh size (h is then the fluid's); empty = h.
    % A closed Sigma is then sampled at each side's own spacing, so the
    % two discrete interfaces differ.
    options.h_solid double {mustBeScalarOrEmpty, mustBePositive} = []
    % Unfitted families only (ignored for fitted). Passed to split only
    % when set; empty leaves split's own policy (45 for interior, 90 = off
    % for the protrude families).
    options.max_normal_deviation double {mustBeScalarOrEmpty} = []
    % Unfitted families only: the same bound for the FLUID side alone;
    % empty (the default) keeps max_normal_deviation on both sides.
    options.max_normal_deviation_fluid double {mustBeScalarOrEmpty} = []
    % Unfitted families only: split's prune_hanging for the FLUID side
    % alone (false keeps the fluid's "teeth"); empty = split's own policy.
    options.prune_hanging_fluid logical {mustBeScalarOrEmpty} = []
    % Unfitted families only. Empty: the crossed grid of step h. A handle:
    % mesh2d meshes `box` to this size target instead (h is then only the
    % boundary default of the mesher), so the background can be graded
    % towards Sigma while remaining unfitted to it.
    options.size_function function_handle = function_handle.empty
end
box = options.box;
closed = is_closed(curve);
fluid_first = ismember(options.fluid_side, ["outer", "upper"]);   % split's mesh_1

switch family
    case "fitted"
        % (h_1, h_2): the sizes of split's (mesh_1, mesh_2)
        h_solid = options.h_solid;  if isempty(h_solid), h_solid = h; end
        h_12 = [h, h_solid];  if ~fluid_first, h_12 = [h_solid, h]; end
        if closed
            [m1, m2] = fitted_closed(curve, box, h_12(1), options.gap, h_12(2));
        else
            assert(isequal(box, [0 1 0 1]), "interface_meshes:box", ...
                "fitted_split meshes graph curves on the unit square only.");
            [m1, m2] = fitted_split(curve, h_12, options.gap);
        end
    otherwise
        % cut elements: to the side that protrudes
        mode = "discard";
        if family ~= "interior"
            to_fluid = family == "protrude_fluid";
            mode = ternary(to_fluid == fluid_first, "mesh_1", "mesh_2");
        end
        if isempty(options.size_function)
            Nx = round((box(2) - box(1)) / h);  Ny = round((box(4) - box(3)) / h);
            bg = StructuredSimplicialGrid(Nx, Ny, box(1), box(2), box(3), box(4), ...
                style = "crossed");
        else
            corners = [box(1) box(2) box(2) box(1); box(3) box(3) box(4) box(4)];
            bg = SimplicialGrid(corners, [], mesh_generator = "mesh2d", ...
                element_diameter = h, size_function = options.size_function);
        end
        split_args = {};
        if ~isempty(options.max_normal_deviation)
            split_args = {"max_normal_deviation", options.max_normal_deviation};
        end
        if ~isempty(options.max_normal_deviation_fluid)
            assert(fluid_first, "interface_meshes:fluidBound", ...
                "max_normal_deviation_fluid needs the fluid on split's mesh_1 (fluid_side outer/upper).");
            split_args = [split_args, {"max_normal_deviation_mesh_1", options.max_normal_deviation_fluid}];
        end
        if ~isempty(options.prune_hanging_fluid)
            assert(fluid_first, "interface_meshes:fluidBound", ...
                "prune_hanging_fluid needs the fluid on split's mesh_1 (fluid_side outer/upper).");
            split_args = [split_args, {"prune_hanging_mesh_1", options.prune_hanging_fluid}];
        end
        [m1, m2] = bg.split(curve, "cut_elements", mode, split_args{:});
end
if fluid_first
    fluid_mesh = m1;  solid_mesh = m2;
else
    fluid_mesh = m2;  solid_mesh = m1;
end
if nargout > 2
    info = struct("family", family, "h", h, ...
        "n_shared", numel(intersect(face_set(fluid_mesh), face_set(solid_mesh))), ...
        "area_gap", (box(2) - box(1))*(box(4) - box(3)) ...
                    - sum(fluid_mesh.volume) - sum(solid_mesh.volume));
end
end

%=========================================================================%
function tf = is_closed(curve)
tf = isa(curve, "LevelSetCurve") || ...
     (isprop(curve, "is_periodic") && curve.is_periodic) || ...
     norm(curve.get_2d_points(curve.t_range(1)) - curve.get_2d_points(curve.t_range(2))) < 1e-10;
end

function [outer_mesh, inner_mesh] = fitted_closed(curve, box, h, gap, h_inner)
% Both sides meshed to the linear interpolant of Sigma: fitted_split's
% closed construction, generalized from the unit square to any box.
% h_inner (default h) sizes the inner side; when it differs from h, the
% inner polygon is the interpolant at its own spacing.
if nargin < 5, h_inner = h; end
P = closed_interpolant(curve, h);  N = size(P, 2);
P_in = P;
if h_inner ~= h
    P_in = closed_interpolant(curve, h_inner);
end
N_in = size(P_in, 2);
% an explicit gap retreats each polygon along the outward normal
if gap > 0
    P_outer = P + gap*polygon_normals(P);  P_inner = P_in - gap*polygon_normals(P_in);
else
    P_outer = P;  P_inner = P_in;
end
corners = [box(1) box(2) box(2) box(1); box(3) box(3) box(4) box(4)];
iface_conn = [4 + [(1:N); (2:N), 1]; 3*ones(1, N)];       % marker 3: Interface
outer_conn = [[1; 2; 1], [2; 3; 1], [3; 4; 1], [4; 1; 1], iface_conn];
outer_mesh = SimplicialGrid([corners, P_outer], outer_conn, ...
    element_diameter = h, mesh_generator = "mesh2d");
inner_conn = [[(1:N_in); (2:N_in), 1]; 3*ones(1, N_in)];
% mesh2d may still insert points on a thin polygon's boundary, so a
% fitted pair is zero-gap but not always face-conforming; the induced
% partition absorbs that and `info.n_shared` reports it.
% refine_boundaries = false is not the remedy: it switches off
% refinement altogether.
inner_mesh = SimplicialGrid(P_inner, inner_conn, ...
    element_diameter = h_inner, mesh_generator = "mesh2d", refine_boundaries = true);
end

function P = closed_interpolant(curve, h)
% Sampling is uniform in arc length, so parametrization speed does not
% cluster the vertices; a LevelSetCurve's induced parametrization is arc
% length, so its uniform-t grid is kept. t_range(1) stays the first
% station, anchoring any distinguished point (an airfoil's tip) there.
% closed_box_split is a separate family (exact-normal, one-sided gaps).
tt_dense = linspace(curve.t_range(1), curve.t_range(2), 4001);
P_dense = curve.get_2d_points(tt_dense);
seg_len = vecnorm(diff(P_dense, 1, 2));
perimeter = sum(seg_len);
N = max(16, round(perimeter / h));
if isa(curve, "LevelSetCurve")
    tt = linspace(curve.t_range(1), curve.t_range(2), N + 1);
else
    s_dense = [0, cumsum(seg_len)];
    tt = interp1(s_dense, tt_dense, linspace(0, perimeter, N + 1));
end
tt = tt(1:end-1);
P = curve.get_2d_points(tt);
% minimum-spacing floor, so near-duplicate vertices cannot reach mesh2d
keep = true(1, size(P, 2));  last = 1;
for j = 2:size(P, 2)
    if norm(P(:, j) - P(:, last)) < 0.2*h, keep(j) = false; else, last = j; end
end
if norm(P(:, find(keep, 1, "last")) - P(:, 1)) < 0.2*h && nnz(keep) > 3
    keep(find(keep, 1, "last")) = false;
end
P = P(:, keep);
end

function n_out = polygon_normals(P)
% the outward unit normal at each vertex of the closed polygon P (of the
% edge that ends there)
n_out = P - circshift(P, 1, 2);  n_out = [n_out(2, :); -n_out(1, :)];
n_out = n_out ./ vecnorm(n_out);
% clockwise polygon (negative signed area): flip to the outward side
if sum(P(1, :).*circshift(P(2, :), -1) - circshift(P(1, :), -1).*P(2, :)) < 0
    n_out = -n_out;
end
end

function s = face_set(m)
% one key per edge, from its two endpoints in a canonical order (the
% lexicographically smaller first) -- a componentwise min/max would be
% the edge's bounding box and equate (0,1)-(1,0) with (0,0)-(1,1)
E = m.faces(1:2, m.intfaces(:));
% rounded to the key's own resolution first, so the ordering decision and
% the key agree on which points are equal
P1 = round(m.coordinates(:, E(1, :)).', 9);  P2 = round(m.coordinates(:, E(2, :)).', 9);
swap = P1(:, 1) > P2(:, 1) | (P1(:, 1) == P2(:, 1) & P1(:, 2) > P2(:, 2));
A = P1;  B = P2;  A(swap, :) = P2(swap, :);  B(swap, :) = P1(swap, :);
s = string(compose("%.9f_%.9f_%.9f_%.9f", [A, B]));
end

function v = ternary(c, a, b)
if c, v = a; else, v = b; end
end
