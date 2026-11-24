function [outer_mesh, inner_mesh] = closed_box_split(interface, h, gap, options)
%-------------------------------------------------------------------------
%Description: Splits an axis-aligned box into two subdomains meeting at a
%             closed interface, meshing each as its own polygon so that its
%             boundary follows the linear interpolation of Sigma.
%
%             The two returned meshes are always ordered
%                 outer_mesh  the region between Sigma and the box
%                 inner_mesh  the region enclosed by Sigma.
%
%             Differs from fitted_split, which does the same job for the
%             unit square: the box is arbitrary, the interface may be any
%             AbstractCurve (including a LevelSetCurve), and polygon points
%             are spaced uniformly in arc length, perimeter/h of them
%             rather than 1/h. The two point counts define different mesh
%             families, so the functions are kept separate.
%
%             Also differs from interface_meshes' fitted family on a
%             closed curve, which samples the same arc-length polygon but
%             opens a symmetric gap along the polygon's edge normals after
%             the spacing floor. Here the gap follows the curve's exact
%             normal, before the floor, and may be one-sided; the two
%             constructions are distinct mesh families.
%
%Usage:       [outer, inner] = closed_box_split(get_interface("kidney"), ...
%                                  1/12, 0, box = [-1 1 -1 1])
%
%See also: fitted_split, generate_mesh_list, SimplicialGrid.split
%-------------------------------------------------------------------------
arguments
    % Any closed curve: both the ParametrizedCurve and the LevelSetCurve
    % lineages are reached through the AbstractCurve query surface only.
    interface AbstractCurve

    % Mesh size. A scalar applies to both sides; [h_outer, h_inner] sizes
    % them independently.
    h double = .1

    % Gap, along the outward normal. A scalar offsets both polygons by that
    % much, leaving Sigma between them at total separation 2*gap;
    % [d_outer, d_inner] offsets them independently, so one side can retreat
    % while the other stays on Sigma. Valid only while the offset stays
    % below the local radius of curvature.
    gap double = 0

    options.box (1,4) double = [0 1 0 1]
    options.boundary_markers double = [1, 1]

    % Polygon points around the whole interface, as a function of h. Empty
    % takes max(16, perimeter/h) -- the 16 keeps a coarse level from
    % degenerating into a triangle.
    options.interface_points function_handle = function_handle.empty

    % Drop a polygon point closer than this fraction of h to its kept
    % predecessor. Without it, a rounded trailing edge hands mesh2d
    % near-duplicate vertices.
    options.min_spacing (1,1) double = 0.2
end
if isscalar(h),   h   = [h, h];     end
if isscalar(gap), gap = [gap, gap]; end
assert(is_periodic(interface), "closed_box_split:openInterface", ...
    "The interface must be closed; use fitted_split for a graph.");

box = options.box;
corners = [box(1) box(2) box(2) box(1); box(3) box(3) box(4) box(4)];
bdy_marker = options.boundary_markers(1);

%% Outer: the box, with the offset interpolant as an internal hole --------
P_outer = offset_polygon(interface, h(1), gap(1), options);
N = size(P_outer, 2);
hole_conn = [4 + [(1:N); (2:N), 1]; 3*ones(1, N)];       % marker 3: Interface
outer_conn = [[1; 2; bdy_marker], [2; 3; bdy_marker], ...
              [3; 4; bdy_marker], [4; 1; bdy_marker], hole_conn];
outer_mesh = SimplicialGrid([corners, P_outer], outer_conn, ...
    element_diameter = h(1), mesh_generator = "mesh2d");

%% Inner: the offset interpolant alone ------------------------------------
P_inner = offset_polygon(interface, h(2), -gap(2), options);
N = size(P_inner, 2);
inner_conn = [[(1:N); (2:N), 1]; 3*ones(1, N)];
inner_mesh = SimplicialGrid(P_inner, inner_conn, ...
    element_diameter = h(2), mesh_generator = "mesh2d", ...
    refine_boundaries = true);
end

%=========================================================================%
function P = offset_polygon(I, h, offset, options)
% Vertices of the interpolating polygon, displaced by `offset` along the
% outward normal (negative displaces inward).
tt_dense = linspace(I.t_range(1), I.t_range(2), 2001);
P_dense = I.get_2d_points(tt_dense);
seg_len = vecnorm(diff(P_dense, 1, 2));
perimeter = sum(seg_len);

if isempty(options.interface_points)
    N = max(16, round(perimeter / h));
else
    N = options.interface_points(h);
end

if isa(I, "LevelSetCurve")
    % Its induced parametrization is already arc length; the uniform-t
    % grid also keeps the seam anchor t_range(1) a vertex.
    tt = linspace(I.t_range(1), I.t_range(2), N + 1);
else
    s_dense = [0, cumsum(seg_len)];
    tt = interp1(s_dense, tt_dense, linspace(0, perimeter, N + 1));
end
tt = tt(1:end-1);                       % the closing point repeats the first

P = I.get_2d_points(tt);
if offset ~= 0
    % The rotated-tangent normal points inward on a clockwise
    % parametrization: orient it against the polygon's own outward normal
    % (central chords, flipped by the signed area) before offsetting.
    n = I.unit_normal(tt);
    d = circshift(P, -1, 2) - circshift(P, 1, 2);
    n_poly = [d(2, :); -d(1, :)];
    if sum(P(1, :).*circshift(P(2, :), -1) - circshift(P(1, :), -1).*P(2, :)) < 0
        n_poly = -n_poly;
    end
    if sum(dot(n, n_poly)) < 0
        n = -n;
    end
    P = P + offset * n;
end

% Minimum-spacing floor, cyclic: drop any point too close to its kept
% predecessor, then re-check the wrap-around pair.
keep = true(1, size(P, 2));
last = 1;
for j = 2:size(P, 2)
    if norm(P(:, j) - P(:, last)) < options.min_spacing * h
        keep(j) = false;
    else
        last = j;
    end
end
if nnz(keep) > 3 && ...
        norm(P(:, find(keep, 1, "last")) - P(:, 1)) < options.min_spacing * h
    keep(find(keep, 1, "last")) = false;
end
P = P(:, keep);
end
