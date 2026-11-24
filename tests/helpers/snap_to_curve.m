function [mesh, info] = snap_to_curve(mesh, curve, options)
%SNAP_TO_CURVE  Interface vertices of a carved mesh moved just outside a curve.
%
%   [mesh, info] = snap_to_curve(mesh, curve, normal_sign = 1)
%
%   The curve's unit normal must point into the mesh's side (the fluid:
%   upward over a graph, outward around a closed curve); normal_sign = -1
%   uses the opposite normal, for the mesh on the other side (the solid:
%   below a graph, inside a closed curve). Each vertex of
%   the faces tagged "Interface" goes to its closest curve point x(t_v),
%   offset along the normal n(t_v) by delta_v, the largest bulge of the
%   curve beyond the chords of its adjacent interface faces (sampled
%   between their end parameters), so that every interface chord lies on
%   the mesh's side of the curve, touching it where it bulges most: the
%   mesh lies inside its exact domain and its transfer paths are
%   O(h^2). Closest-point parameters that run against the order of the
%   interface chain (vertices near a corner of the curve projected past
%   each other, which would fold an element) are replaced by evenly spaced
%   ones between the nearest in-order neighbours. The generalization of
%   snap_outside (circles, radial moves) to any curve with a closest-point
%   map and a normal. Vertices mapped to
%   an end of an open curve are put on that end, with no offset (on a
%   WallEndedCurve they stay on the wall). A vertex whose move would
%   invert an element is left unsnapped, with the other interface
%   vertices of that element.
%
%   info (optional): unsnapped, the interface vertices (mesh indices) left
%   where the carve put them, and n_unsnapped, their number; n_reordered,
%   the interface vertices whose closest-point parameter was replaced by
%   the chain ordering; n_unordered, the interface faces whose end
%   parameters still run against the chain (no in-order neighbours within
%   the search window); max_distance, the largest distance of an
%   interface vertex of the returned mesh to the curve (closest point).
%
%See also: snap_outside, lagrangian_carved_march, WallEndedCurve
arguments
    mesh
    curve
    options.normal_sign (1,1) double {mustBeMember(options.normal_sign, [-1 1])} = 1
end
X = mesh.coordinates;
F = mesh.faces(1:2, mesh.intfaces);
[v, ~, e] = unique(F(:));
e = reshape(e, 2, []);
t = reshape(curve.map_points_to_interface(X(:, v)), 1, []);
period = diff(curve.t_range);
[t, n_reordered, n_unordered] = order_parameters(mesh, v, e, t, curve.t_range, is_periodic(curve));
P = curve.get_2d_points(t);
N = options.normal_sign*reshape(curve.unit_normal(t), 2, []);
delta = zeros(1, numel(v));
theta = linspace(0, 1, 11);
for f = 1:size(e, 2)
    a = e(1, f);  b = e(2, f);
    tb = t(b);
    if is_periodic(curve) && abs(tb - t(a)) > period/2
        tb = tb - sign(tb - t(a))*period;   % the short way round
    end
    C = curve.get_2d_points(t(a) + theta*(tb - t(a)));
    chord = P(:, a) + theta .* (P(:, b) - P(:, a));
    n = (N(:, a) + N(:, b))/2;
    bulge = max(sum((C - chord) .* n, 1));
    delta([a, b]) = max(delta([a, b]), bulge);
end
if ~is_periodic(curve)
    delta(t == curve.t_range(1) | t == curve.t_range(2)) = 0;
end
X0 = X;
X(:, v) = P + delta .* N;
% a vertex whose move would invert an element (a carve that left a gap
% of more than about one element, e.g. after the hanging-element prune)
% stays where the carve put it, as do the other interface vertices of
% that element; repeated until no element is inverted
area = @(X) .5*sum(reshape(X(1, mesh.elements), 3, []).*circshift(reshape(X(2, mesh.elements), 3, []), -1) ...
    - circshift(reshape(X(1, mesh.elements), 3, []), -1).*reshape(X(2, mesh.elements), 3, []));
inverted = area(X) <= 0;
while any(inverted)
    back = intersect(unique(mesh.elements(:, inverted)), v);
    X(:, back) = X0(:, back);
    inverted = area(X) <= 0;
end
mesh.coordinates = X;
mesh.volume = area(X);
dx_e = diff(reshape(X(1, mesh.faces(1:2, :)), 2, []));
dy_e = diff(reshape(X(2, mesh.faces(1:2, :)), 2, []));
mesh.area = sqrt(dx_e.^2 + dy_e.^2);
n_by_face = [dy_e; -dx_e];
mesh.normals = reshape(n_by_face(:, mesh.boundaries), 2*3, []);
if nargout > 1
    unsnapped = v(all(X(:, v) == X0(:, v), 1) & any(X0(:, v) ~= P + delta .* N, 1));
    t_final = reshape(curve.map_points_to_interface(X(:, v)), 1, []);
    info = struct("unsnapped", reshape(unsnapped, 1, []), "n_unsnapped", numel(unsnapped), ...
        "n_reordered", n_reordered, "n_unordered", n_unordered, ...
        "max_distance", max(vecnorm(X(:, v) - curve.get_2d_points(t_final)), [], "all"));
end
end

function [t, n_reordered, n_unordered] = order_parameters(mesh, v, e, t, t_range, periodic)
% the closest-point parameters t of the interface vertices v, made to run
% along the curve in the order of the interface chain: the faces e (local
% indices into v) are oriented with the mesh on their left, and a face
% whose end parameters run against the chain's prevailing direction (two
% vertices near a corner of the curve projected past each other) gets its
% two ends, and if needed more chain neighbours on both sides, replaced
% by parameters evenly spaced between the nearest pair of vertices before
% and after it that are in order. n_reordered counts the vertices whose
% parameter was replaced, n_unordered the faces left against the chain
max_window = 8;
period = diff(t_range);
T = mesh.elements;
n = numel(v);
global_e = v(e);
for f = 1:size(e, 2)                     % orient each face by its element
    k = find(any(T == global_e(1, f), 1) & any(T == global_e(2, f), 1), 1);
    c = setdiff(T(:, k), global_e(:, f));
    A = mesh.coordinates(:, [global_e(:, f); c]);
    if det([A(:, 2) - A(:, 1), A(:, 3) - A(:, 1)]) < 0
        e(:, f) = e([2 1], f);
    end
end
next = zeros(1, n);  prev = zeros(1, n);
next(e(1, :)) = e(2, :);  prev(e(2, :)) = e(1, :);
% (step reads t when called: the windows below change it)
step = @(t, a, b) wrap(t(b) - t(a), period, periodic);
direction = sign(sum(step(t, e(1, :), e(2, :))));
t0 = t;
for f = find(direction*step(t, e(1, :), e(2, :)) <= 0)
    p = e(1, f);  q = e(2, f);
    if direction*step(t, p, q) > 0, continue, end   % fixed by an earlier window
    for grown = 1:max_window
        if direction*step(t, p, q) > 0 || prev(p) == 0 || next(q) == 0, break, end
        p = prev(p);  q = next(q);
    end
    if direction*step(t, p, q) <= 0, continue, end
    chain = p;
    while chain(end) ~= q
        chain(end + 1) = next(chain(end)); %#ok<AGROW>
    end
    t(chain(2:end - 1)) = t(p) + step(t, p, q)*(1:numel(chain) - 2)/(numel(chain) - 1);
    if periodic
        t(chain) = t_range(1) + mod(t(chain) - t_range(1), period);
    end
end
n_reordered = nnz(t ~= t0);
n_unordered = nnz(direction*wrap(t(e(2, :)) - t(e(1, :)), period, periodic) <= 0);
end

function d = wrap(d, period, periodic)
% a parameter difference the short way round a periodic curve
if periodic
    d = d - period*round(d/period);
end
end
