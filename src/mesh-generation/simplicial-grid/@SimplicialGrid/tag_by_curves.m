function [obj, assignment] = tag_by_curves(obj, curves, tags, options)
%-------------------------------------------------------------------------
%Description: Assigns every boundary face to one of the curves C_1..C_N by
%             closest-point projection of its quadrature points, and tags
%             it with that curve's tag. The curves must partition the
%             boundary (disjoint except at endpoints): a curve that
%             contains another (a full circle and an arc of it) makes the
%             faces on the shared stretch ties.
%
%             A face is resolved when every quadrature point has the same
%             nearest curve and the runner-up trails it by at least
%             margin*|F| at every point. Otherwise it is unresolved: its
%             curve is NaN and it keeps its old tag. A resolved face is
%             still unresolved when
%               "far"   a point is farther than reach*|F| from every curve
%               "side"  a point's closest curve point lies behind the face
%                       (against its outward normal, beyond 60 degrees),
%                       farther than behind*|F|: the projection has jumped
%                       to the wrong stretch of a thin geometry.
%             Faces whose projections tie at a junction of two curves are
%             unresolved by construction ("ambiguous"). Remedies, after
%             which the call is repeated: split the face where it crosses
%             the junction (bisect_boundary_face), or assign it by hand
%             (force).
%
%             Comparing with the existing tagging: assignment.previous_tag
%             holds faces(3,:) on entry and assignment.changed lists the
%             resolved faces whose tag changed.
%
%Usage:       [mesh, a] = tag_by_curves(mesh, {flag, root}, ["Interface", "Clamp"])
%             unresolved = find(isnan(a.curve));
%             mesh = tag_by_curves(mesh, curves, tags, force = [faces; curve_ids])
%
%             assignment fields (1 x n_faces unless noted):
%               curve         curve index; 0 = not judged (interior or
%                             not in options.faces); NaN = unresolved
%               status        "" | "ok" | "majority" | "forced" |
%                             "ambiguous" | "far" | "side"
%               distance      2 x n_faces: max over the points of the
%                             nearest distance; min of the runner-up's
%               margin        min over the points of (runner-up - nearest),
%                             in units of |F|
%               candidates    2 x n_faces: the two curves with the smallest
%                             mean distance
%               previous_tag  faces(3,:) on entry
%               changed       resolved faces whose tag changed
%
%See also: bisect_boundary_face, register_tag, AbstractCurve.distance
%-------------------------------------------------------------------------
arguments
    obj SimplicialGrid
    curves (1, :) cell                    % AbstractCurve each
    tags (1, :) string                    % tag per curve (may repeat)
    options.faces (1, :) double = []      % faces to judge; empty = all boundary faces
    options.order (1,1) double = 2        % order + 1 Gauss points per face
    options.margin (1,1) double = 0.05
    options.reach (1,1) double = 5
    options.behind (1,1) double = 0.25   % Inf disables the "side" check
    options.force (2, :) double = zeros(2, 0)   % [face; curve index], applied last
    options.conflicts (1,1) string {mustBeMember(options.conflicts, ["mark", "majority"])} = "mark"
end
n_curves = numel(curves);
assert(numel(tags) == n_curves, "SimplicialGrid:tag_by_curves:tags", ...
    "one tag per curve: %d curves, %d tags", n_curves, numel(tags));
assert(all(cellfun(@(c) isa(c, "AbstractCurve"), curves)), ...
    "SimplicialGrid:tag_by_curves:curves", "every curve must be an AbstractCurve");
n_faces = size(obj.faces, 2);

% boundary faces: incidence lists them twice with the same (element, side)
K = permute(obj.incidence(:, 1, :), [1 3 2]);
s = permute(obj.incidence(:, 2, :), [1 3 2]);
is_boundary = K(1, :) == K(2, :) & s(1, :) == s(2, :);
faces = options.faces;
if isempty(faces)
    faces = find(is_boundary);
end
assert(all(is_boundary(faces)), "SimplicialGrid:tag_by_curves:interior", ...
    "only boundary faces can be tagged");

warn_overlap(curves);

% quadrature points (2 x nq x nf) and outward normals of the judged faces
[s_q, ~] = get_1d_reference_quadrature(options.order);
nq = numel(s_q);
nf = numel(faces);
X = get_1d_physical(obj, repmat(s_q(:), 1, nf), faces);
X = reshape(X, 2, []);
N = obj.exterior_unit_normals;
N = reshape(N(:, 1, sub2ind([3, size(obj.elements, 2)], s(1, faces), K(1, faces))), 2, nf);
F = obj.area(faces);

% distances to every curve: D(j, q, f), closest points P(:, q, f, j)
D = zeros(n_curves, nq*nf);
P = zeros(2, nq*nf, n_curves);
for j = 1:n_curves
    [D(j, :), ~, P(:, :, j)] = distance(curves{j}, X);
end
[Ds, order] = sort(D, 1);
d1 = reshape(Ds(1, :), nq, nf);
if n_curves > 1
    d2 = reshape(Ds(2, :), nq, nf);
else
    d2 = Inf(nq, nf);
end
nearest = reshape(order(1, :), nq, nf);

curve = nan(1, nf);
status = repmat("ambiguous", 1, nf);
agree = all(nearest == nearest(1, :), 1);
gap = min(d2 - d1, [], 1) ./ F;
clean = agree & gap >= options.margin;
curve(clean) = nearest(1, clean);
status(clean) = "ok";
if options.conflicts == "majority"
    [mode_curve, count] = mode(nearest, 1);
    win = ~clean & count >= 2/3*nq;
    curve(win) = mode_curve(win);
    status(win) = "majority";
end

% far: beyond the reach of every curve
far = ~isnan(curve) & any(d1 > options.reach*F, 1);
% side: the closest point of the assigned curve lies behind the face
curve_pt = reshape(repmat(curve, nq, 1), 1, []);
P_own = X;
for j = 1:n_curves
    on_j = curve_pt == j;
    P_own(:, on_j) = P(:, on_j, j);
end
V = reshape(P_own - X, 2, nq, nf);
d_own = reshape(vecnorm(reshape(V, 2, [])), nq, nf);
cos_angle = reshape(sum(V .* reshape(N, 2, 1, nf), 1), nq, nf) ./ max(d_own, realmin);
side = ~isnan(curve) & any(cos_angle < -cosd(60) & d_own > options.behind*F, 1);
curve(far) = NaN;   status(far) = "far";
curve(side & ~far) = NaN;   status(side & ~far) = "side";

% into full-length arrays
assignment.curve = zeros(1, n_faces);
assignment.curve(faces) = curve;
assignment.status = repmat("", 1, n_faces);
assignment.status(faces) = status;
assignment.distance = nan(2, n_faces);
assignment.distance(:, faces) = [max(d1, [], 1); min(d2, [], 1)];
assignment.margin = nan(1, n_faces);
assignment.margin(faces) = gap;
mean_D = reshape(mean(reshape(D, n_curves, nq, nf), 2), n_curves, nf);
[~, rank] = sort(mean_D, 1);
assignment.candidates = zeros(2, n_faces);
assignment.candidates(1:min(2, n_curves), faces) = rank(1:min(2, n_curves), :);

% manual assignment, applied last
forced = options.force(1, :);
assert(all(is_boundary(forced)), "SimplicialGrid:tag_by_curves:interior", ...
    "only boundary faces can be forced");
assert(all(ismember(options.force(2, :), 1:n_curves)), ...
    "SimplicialGrid:tag_by_curves:forceCurve", "forced curve index out of range");
assignment.curve(forced) = options.force(2, :);
assignment.status(forced) = "forced";

% tags of the resolved faces
assignment.previous_tag = obj.faces(3, :);
tag_number = zeros(1, n_curves);
for j = 1:n_curves
    [obj, tag_number(j)] = register_tag(obj, tags(j));
end
resolved = find(assignment.curve > 0);
obj.faces(3, resolved) = tag_number(assignment.curve(resolved));
assignment.changed = resolved(obj.faces(3, resolved) ~= assignment.previous_tag(resolved));
end

%=========================================================================%
function warn_overlap(curves)
% Curves sharing more than their endpoints: samples of one lying on another
n_samples = 200;
n = numel(curves);
S = cell(1, n);
for j = 1:n
    t = linspace(curves{j}.t_range(1), curves{j}.t_range(2), n_samples);
    S{j} = reshape(get_2d_points(curves{j}, t), 2, []);
end
for a = 1:n
    for b = [1:a-1, a+1:n]
        scale = max(vecnorm(S{b} - S{b}(:, 1)));
        spacing = max(vecnorm(diff(S{b}, 1, 2)));
        % only samples near b's sampling are projected
        near = min(vecnorm(S{a} - permute(S{b}, [1 3 2])), [], 3) <= spacing;
        if nnz(near) <= 2, continue; end
        d = distance(curves{b}, S{a}(:, near));
        if nnz(d <= 1e-6*max(scale, eps)) > 2
            warning("SimplicialGrid:tag_by_curves:overlappingCurves", ...
                "curves %d and %d overlap on more than their endpoints", a, b);
        end
    end
end
end
