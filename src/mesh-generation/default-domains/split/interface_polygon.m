function poly = interface_polygon(mesh, curve)
%INTERFACE_POLYGON  A mesh's own interface facets, as a PolygonalCurve.
%
%   poly = interface_polygon(mesh, curve)
%
%   The interface mapping of a naive comparison: hand the solver the
%   mesh's interface polyline as if it were the true Sigma, so every
%   transfer path is zero-length and whatever the facets get wrong about
%   the geometry enters the solution directly. Closed if the facets form
%   a loop, open otherwise. `curve` fixes the orientation to the true
%   Sigma's (closed: same signed area; open: same direction of travel),
%   so the normal the coupling builds from it points the same way.
%
%See also: interface_meshes, PolygonalCurve
arguments
    mesh  SimplicialGrid
    curve AbstractCurve
end
E = mesh.faces(1:2, mesh.intfaces(:));
n_e = size(E, 2);
[vid, ~, ic] = unique(E(:));
deg = accumarray(ic, 1);
ends = vid(deg == 1);
closed = isempty(ends);
assert(closed || numel(ends) == 2, "interface_polygon:notAPolyline", ...
    "the interface facets form neither a loop nor a single chain (%d endpoints).", numel(ends));
% walk
if closed, v = E(1, 1); else, v = ends(1); end
order = zeros(1, n_e + ~closed);  used = false(1, n_e);
for i = 1:n_e
    order(i) = v;
    e = find(~used & any(E == v, 1), 1);
    used(e) = true;
    v = E(E(:, e) ~= v, e);
end
if ~closed, order(end) = v; end
V = mesh.coordinates(:, order);
% orientation
tt = linspace(curve.t_range(1), curve.t_range(2), 64);
ref = curve.get_2d_points(tt);
if closed
    sa = @(P) 0.5*sum(P(1, :).*circshift(P(2, :), -1) - circshift(P(1, :), -1).*P(2, :));
    if sign(sa(V)) ~= sign(sa(ref(:, 1:end-1))), V = fliplr(V); end
else
    if dot(V(:, end) - V(:, 1), ref(:, end) - ref(:, 1)) < 0, V = fliplr(V); end
end
% normal orientation: agree with the curve's unit_normal, compared at the
% segment midpoints against the nearest curve sample
poly = PolygonalCurve(V, closed = closed);
t_mid = (0:size(poly.chords, 2) - 1) + 0.5;
mid = poly.get_2d_points(t_mid);
[~, nearest] = min((mid(1, :).' - ref(1, :)).^2 + (mid(2, :).' - ref(2, :)).^2, [], 2);
agreement = sum(poly.unit_normal(t_mid) .* curve.unit_normal(tt(nearest)), "all");
if agreement < 0
    poly = PolygonalCurve(V, closed = closed, normal_sign = -1);
end
end
