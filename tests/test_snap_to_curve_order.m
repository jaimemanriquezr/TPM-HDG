% test_snap_to_curve_order.m
% snap_to_curve on a carved boundary that overhangs (a pruned gap next to
% a corner of the curve): the interface chain 1-2-3-4-5 runs back to the
% left between vertices 2 and 3, so their closest points on y = 0.5 are
% out of chain order. Snapping them there would fold element (2, 3, 9);
% the chain-ordered parameters put every interface vertex on the curve,
% in order, with no element inverted:
%   1. all interface vertices on y = 0.5, x increasing along the chain;
%   2. vertices 2 and 3 evenly spaced between their in-order neighbours;
%   3. every element keeps a positive area;
%   4. an ordinary carve of a circle: only the interface vertices move,
%      onto the curve or just outside it (offset below h^2);
%   5. the second output, info: no vertex unsnapped in 1-4 (two reordered
%      in the chain), and a curve above the chain's interior vertices
%      forces the fallback, reported vertex by vertex.
%   6. normal_sign picks the side: on the circle of check 4 and on the graph
%      y = 0.5 + 0.05 sin(2 pi x), the interface vertices of the fluid carve
%      (outside, above; normal_sign = 1) lie on or outside the curve and those
%      of the solid carve (inside, below; normal_sign = -1) on or inside it
%      (round-off 1e-14); the solid carve snapped with normal_sign = 1 does
%      not (on a flat curve both signs give the same mesh).
clc; clear;

n_fail = 0;

%% 1-3. the overhanging chain
X = [0 .55; .4 .56; .3 .66; .7 .6; 1 .55; 0 1; .5 1; 1 1; .15 .8; .6 .8].';
T = [1 2 9; 2 3 9; 1 9 6; 9 3 7; 9 7 6; 3 4 10; 3 10 7; 4 5 10; 10 5 8; 10 8 7].';
mesh = SimplicialGrid(X, T, mesh_generator = "array");
mesh = apply_tag(mesh, "Interface", @(P) (P(2, :) < 0.7 & P(1, :) > 1e-12 & P(1, :) < 1 - 1e-12) ...
    | (P(2, :) < 0.62 & (P(1, :) < 1e-12 | P(1, :) > 1 - 1e-12)));
n_fail = n_fail + check(isequal(sort(mesh.faces(1:2, mesh.intfaces), 1), [1 2 3 4; 2 3 4 5]), ...
    "the interface chain is not 1-2-3-4-5");
snapped = snap_to_curve(mesh, FlatCurveHorizontal(0.5));
Y = snapped.coordinates(:, 1:5);
fprintf("snapped chain x: %s; y - 0.5 up to %.1e; min area %.2e\n", ...
    sprintf("%.4f ", Y(1, :)), max(abs(Y(2, :) - 0.5)), min(snapped.volume));
n_fail = n_fail + check(max(abs(Y(2, :) - 0.5)) < 1e-14, "an interface vertex is off the curve");
n_fail = n_fail + check(all(diff(Y(1, :)) > 0), "the snapped chain is out of order");
n_fail = n_fail + check(max(abs(Y(1, 2:3) - [0.7, 1.4]/3)) < 1e-14, ...
    "vertices 2, 3 not evenly spaced in (0, 0.7)");
n_fail = n_fail + check(all(snapped.volume > 0), "an element is inverted");

%% 4. in-order closest points: only the interface vertices move, onto the circle
circle = CircularCurve([0.5; 0.5], 0.4);
carved = StructuredSimplicialGrid(16, 16, "crossed").split(circle);
snapped = snap_to_curve(carved, circle);
v = unique(carved.faces(1:2, carved.intfaces));
free = setdiff(1:size(carved.coordinates, 2), v);
r = vecnorm(snapped.coordinates(:, v) - [0.5; 0.5]);
fprintf("circle: radius - 0.4 in [%.1e, %.1e]\n", min(r) - 0.4, max(r) - 0.4);
n_fail = n_fail + check(isequal(snapped.coordinates(:, free), carved.coordinates(:, free)), ...
    "circle: a vertex off the interface moved");
n_fail = n_fail + check(min(r) >= 0.4 - 1e-14 && max(r) - 0.4 < 1/16^2, ...
    "circle: interface vertices not just outside the curve");

%% 5. the diagnostics (second output): the cases above, and a forced fallback
[snapped_info, info] = snap_to_curve(carved, circle);
fprintf("circle info: unsnapped %d, reordered %d, unordered %d, max distance %.2e\n", ...
    info.n_unsnapped, info.n_reordered, info.n_unordered, info.max_distance);
n_fail = n_fail + check(isequal(snapped_info.coordinates, snapped.coordinates), ...
    "circle: the second output changed the mesh");
n_fail = n_fail + check(info.n_unsnapped == 0 && info.n_reordered == 0 && info.n_unordered == 0 ...
    && info.max_distance < 1/16^2, "circle: unexpected diagnostics");
[~, info] = snap_to_curve(mesh, FlatCurveHorizontal(0.5));
fprintf("chain info: unsnapped %d, reordered %d, unordered %d, max distance %.1e\n", ...
    info.n_unsnapped, info.n_reordered, info.n_unordered, info.max_distance);
n_fail = n_fail + check(info.n_unsnapped == 0 && info.n_reordered == 2 && info.n_unordered == 0 ...
    && info.max_distance < 1e-14, "chain: unexpected diagnostics");
% y = 0.9 lies above the interior vertices 9 and 10: snapping vertices 1
% and 2 there inverts element (1, 2, 9), and the reverts cascade along the
% chain, so every interface vertex keeps its carve position
[fallback, info] = snap_to_curve(mesh, FlatCurveHorizontal(0.9));
fprintf("fallback info: unsnapped %s, reordered %d, max distance %.4f\n", ...
    mat2str(info.unsnapped), info.n_reordered, info.max_distance);
n_fail = n_fail + check(isequal(fallback.coordinates, snap_to_curve(mesh, FlatCurveHorizontal(0.9)).coordinates), ...
    "fallback: the second output changed the mesh");
n_fail = n_fail + check(isequal(info.unsnapped, 1:5) && info.n_unsnapped == 5 ...
    && isequal(fallback.coordinates(:, 1:5), mesh.coordinates(:, 1:5)), ...
    "fallback: the unsnapped vertices are not reported");
n_fail = n_fail + check(abs(info.max_distance - 0.35) < 1e-14, "fallback: max distance is not 0.9 - 0.55");

%% 6. the side of the snapped vertices: fluid normal_sign = 1, solid -1
grid = StructuredSimplicialGrid(16, 16, "crossed");
graph = @(x) 0.5 + 0.05*sin(2*pi*x);
cases = {"circle", circle, @(X) vecnorm(X - [0.5; 0.5]) - 0.4; ...
         "graph", SmoothCurveHorizontal("0.5 + 0.05*sin(2*pi*t)"), @(X) X(2, :) - graph(X(1, :))};
for c = 1:size(cases, 1)
    [name, curve, side] = cases{c, :};               % side > 0: outside the circle, above the graph
    [fluid, solid] = grid.split(curve);
    interface_side = @(m) side(m.coordinates(:, unique(m.faces(1:2, m.intfaces))));
    s_fluid = interface_side(snap_to_curve(fluid, curve));
    s_solid = interface_side(snap_to_curve(solid, curve, normal_sign = -1));
    s_wrong = interface_side(snap_to_curve(solid, curve, normal_sign = 1));
    fprintf("%s: side of the interface vertices, fluid [%.1e, %.1e], solid [%.1e, %.1e], " + ...
        "solid with normal_sign = 1 [%.1e, %.1e]\n", name, min(s_fluid), max(s_fluid), ...
        min(s_solid), max(s_solid), min(s_wrong), max(s_wrong));
    n_fail = n_fail + check(min(s_fluid) >= -1e-14 && max(s_solid) <= 1e-14, ...
        name + ": snapped interface vertices on the wrong side of the curve");
    n_fail = n_fail + check(max(s_wrong) > 1e-14, name + ": normal_sign does not change the side");
end

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_snap_to_curve_order:failed", "%d check(s) failed", n_fail);
end
