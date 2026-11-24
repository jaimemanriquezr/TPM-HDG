function [s, switched, moved, info] = remap_solid_state(sp, s, old_mesh, new_mesh, k, old_keys, new_keys, ...
    old_face_keys, new_face_keys, options)
%REMAP_SOLID_STATE  Newmark solid history from one carve of a background grid to another.
%
%   [s, switched, moved, info] = remap_solid_state(sp, s, old_mesh, new_mesh, k, ...
%       old_keys, new_keys, old_face_keys, new_face_keys, ...
%       untranslated = [], shift = [0; 0], carry_tolerance = 0)
%
%   s is a Newmark state.solid on old_mesh (LinearElasticityProblem sp,
%   degree k); both meshes are carved from the same background grid, their
%   elements identified by the keys (centroids of the unsnapped carve, one
%   row per element) and their faces by the face keys (midpoints of the
%   unsnapped carve, one row per face). An element present in both carves
%   is the source of itself; an element that switches on takes as source
%   the old element with the nearest centroid. Per field:
%     volume, rotation (Dubiner)  copied on an element of unchanged shape;
%       otherwise the L2 projection of the source's polynomial, evaluated
%       (extrapolated) on the new element, which re-expresses it exactly
%       on a snapped element; the displacement, on an element that
%       switches on or changes shape, is instead the projection of the
%       source's P_k+1 postprocessing u* (postprocessed_displacement): the
%       next step's u_h is close to the projection of u on its element,
%       and a history that differs from it by e (u_h's own polynomial,
%       re-expressed on a new shape or extrapolated, has e of order
%       h^(k+1)) enters the Newmark acceleration as e/(beta dt^2);
%     stress (PEERS)  copied on an element of unchanged shape; otherwise
%       the L2 projection onto PEERS of the new element of the source's
%       stress (exact on P_k tensors: PEERS is not closed under
%       extrapolation, its bubbles depend on the element);
%     skeleton (Legendre)  copied on a face present in both carves with
%       the same end points (submesh keeps the background's face
%       orientation); on a face that is an old face translated by shift
%       (an interface face snapped onto a curve that moved with the
%       solid), the old face's trace before the translation, given as
%       untranslated (the state.solid that translate_state received): the
%       history along the characteristics, with no extension; on any other
%       face the Legendre projection of the trace of the source of an
%       element owning it: for the displacement, its u*, whose trace is one
%       order more accurate than u_h's; for the velocity and acceleration,
%       their volume polynomials, except on a boundary face when
%       untranslated is given: there the old boundary trace slid along the
%       characteristic (slid_traces). The volume trace of the velocity
%       differs from the Newmark trace by O(h^(k+1)) and enters the
%       kinematic row directly; the displacement trace of u* does not need
%       this.
%   slid_subdivisions = m (default 1) replaces the Gauss rule of the
%   rebuilt faces by the same rule on m equal subintervals of each face:
%   the traces read there are only piecewise polynomial along a new face.
%   slid_subdivisions = Inf projects the slid traces exactly instead
%   (slid_traces_exact) where the old boundary polygon is convex near the
%   face: each new boundary face is split where the closest point of the
%   old polygon passes an old vertex, the slid trace being a polynomial on
%   each piece. At a reflex old vertex the closest face also switches
%   where two faces are equally distant, which is not a split point: there
%   the rule is Gauss on a piece with a jump. slid_correction = false drops the
%   volume correction u_K(z) - u_K(p) of the slid traces;
%   slid_displacement = true slides the displacement trace as well.
%   Every field's column count is checked against the new mesh, since a
%   skeleton history of the wrong size is indexed without an error by the
%   interface rows. displacement_previous must be empty (Newmark only).
%   switched = [on, off]: the background elements entering and leaving
%   the solid; moved, the new elements that switch on or change shape.
%
%   carry_tolerance = r > 0 (with untranslated) widens the carried rule
%   to an interface face whose two ends lie within r times its length of
%   the ends, translated by shift, of an old interface face (the nearest
%   one, ends in the same order): an interface snapped onto a curve that
%   moved with the solid only up to a change of its shape. The old trace
%   is then copied unchanged onto a face of slightly different geometry.
%   r = 0 is the exact rule above alone.
%   info (optional, needs untranslated): n_carried, the faces carried by
%   either rule; mismatch, the largest over the new interface faces of
%   the distance of its ends to those of the nearest translated old
%   interface face, over its length (NaN without interface faces).
%
%See also: eulerian_solid_march, LinearElasticityProblem.translate_state
arguments
    sp (1,1) LinearElasticityProblem
    s (1,1) struct
    old_mesh SimplicialGrid
    new_mesh SimplicialGrid
    k (1,1) double
    old_keys (:,2) double
    new_keys (:,2) double
    old_face_keys (:,2) double
    new_face_keys (:,2) double
    options.untranslated = []
    options.shift (2,1) double = [0; 0]
    options.carry_tolerance (1,1) double {mustBeNonnegative} = 0
    % rebuilt faces: Gauss rule on this many equal subintervals of each face
    options.slid_subdivisions (1,1) double {mustBePositive} = 1
    options.slid_correction (1,1) logical = true
    options.slid_displacement (1,1) logical = false
end
assert(options.slid_subdivisions == Inf || mod(options.slid_subdivisions, 1) == 0, ...
    "remap_solid_state:subdivisions", "slid_subdivisions is a positive integer or Inf.");
volume   = ["displacement", "velocity", "acceleration"];
skeleton = volume + "_skeleton";
stress   = ["stress", "stress_velocity", "stress_acceleration"];
rotation = ["rotation", "rotation_velocity", "rotation_acceleration"];
known = [volume, skeleton, stress, rotation, "displacement_previous"];
unknown = setdiff(string(fieldnames(s)).', known);
assert(isempty(unknown), "remap_solid_state:unknownField", "Unknown solid state field(s): %s.", ...
    strjoin(unknown, ", "));
assert(~isfield(s, "displacement_previous") || isempty(s.displacement_previous), ...
    "remap_solid_state:scheme", "displacement_previous is set: Newmark states only.");
has = @(name) isfield(s, name) && isnumeric(s.(name)) && ~isempty(s.(name)) && ~isequal(s.(name), 0);
s_old = s;

% elements: sources and shapes
n_new = size(new_mesh.elements, 2);
[present, from] = ismembertol(new_keys, old_keys, 1e-10, "ByRows", true);
switched = [nnz(~present), size(old_mesh.elements, 2) - nnz(present)];
c_old = centroid(old_mesh);
c_new = centroid(new_mesh);
source = zeros(n_new, 1);
source(present) = from(present);
for K = find(~present).'
    [~, source(K)] = min(vecnorm(c_old - c_new(K, :), 2, 2));
end
kept = present;
kept(present) = all(reshape(abs(new_mesh.coordinates(:, new_mesh.elements(:, present)) ...
    - old_mesh.coordinates(:, old_mesh.elements(:, from(present)))) < 1e-14, 6, []), 1).';
moved = find(~kept);

% volume and rotation fields: Dubiner, component blocks of d2 rows
d2 = nchoosek(k + 2, 2);
% the rule of argument k + 2 is exact to degree 2k + 5: the fits below (onto
% Dubiner P_k and onto PEERS) are exact L2 projections of sources of degree
% k + 1 at most (u_h, u*, the PEERS stress), integrands of degree 2k + 2
[rs, w] = get_2d_reference_quadrature(k + 2);
nq = numel(w);
fit = fit_matrix(dubiner_basis(k, rs), w);
rs_old = zeros(2, nq, 0);
G_new = zeros(2, 2, 0);
G_old = zeros(2, 2, 0);
if ~isempty(moved)
    xy = new_mesh.get_2d_physical(repmat(rs, 1, 1, numel(moved)), moved);
    rs_old = old_mesh.get_2d_local(xy, source(moved));
    G_new = get_affine_map(new_mesh, moved, inverse = true);
    G_old = get_affine_map(old_mesh, source(moved), inverse = true);
end
for name = [volume, rotation]
    if ~has(name), continue, end
    s.(name) = remap_dubiner(s.(name), k, d2, fit, kept, source, moved, rs_old, n_new);
end
% the displacement on an element that switches on or changes shape: the
% projection of the source's u*
U_star = [];
if has("displacement") && ~isempty(moved)
    U_star = postprocessed_displacement(sp, old_mesh, k, s_old);
    d2p = nchoosek(k + 3, 2);
    for j = 1:numel(moved)
        Q = dubiner_basis(k + 1, rs_old(:, :, j));
        for c = 0:1
            s.displacement(c*d2 + (1:d2), moved(j)) = fit*(Q*U_star(c*d2p + (1:d2p), source(moved(j))));
        end
    end
end

% stress fields: PEERS
for name = stress
    if ~has(name), continue, end
    S = s.(name);
    T = zeros(size(S, 1), n_new);
    T(:, kept) = S(:, source(kept));
    for j = 1:numel(moved)
        Phi = peers_points(k, rs, repmat(G_new(:, :, j), 1, 1, nq));           % nq x nb x 4
        Psi = peers_points(k, rs_old(:, :, j), repmat(G_old(:, :, j), 1, 1, nq));
        values = pagemtimes(Psi, S(:, source(moved(j))));                    % nq x 1 x 4
        A = reshape(permute(Phi, [1 3 2]), nq*4, []);                        % (q, c) x nb
        W = repmat(w(:), 4, 1);
        b = reshape(permute(values, [1 3 2]), nq*4, 1);
        % normal equations: the rule's weights are signed, so the problem
        % is not a least-squares one in sqrt(W); solved with the mass
        % matrix scaled to a unit diagonal, since the bubble columns scale
        % with G^2 and the Dubiner ones do not
        M = A.'*(W.*A);
        D = 1./sqrt(diag(M));
        T(:, moved(j)) = D.*((D.*M.*D.') \ (D.*(A.'*(W.*b))));
    end
    s.(name) = T;
end

% skeleton fields: Legendre, component blocks of d1 rows
n_faces = size(new_mesh.faces, 2);
[present_f, from_f] = ismembertol(new_face_keys, old_face_keys, 1e-10, "ByRows", true);
kept_f = present_f;
kept_f(present_f) = all(abs(new_mesh.coordinates(:, new_mesh.faces(1, present_f)) ...
    - old_mesh.coordinates(:, old_mesh.faces(1, from_f(present_f)))) < 1e-14, 1).' ...
    & all(abs(new_mesh.coordinates(:, new_mesh.faces(2, present_f)) ...
    - old_mesh.coordinates(:, old_mesh.faces(2, from_f(present_f)))) < 1e-14, 1).';
carried = false(n_faces, 1);
from_c = zeros(n_faces, 1);
if ~isempty(options.untranslated)
    ends = @(m, F, x) [(m.coordinates(:, m.faces(1, F)) - x).', (m.coordinates(:, m.faces(2, F)) - x).'];
    candidates = find(~kept_f);
    [carried(candidates), from_c(candidates)] = ismembertol(ends(new_mesh, candidates, options.shift), ...
        ends(old_mesh, 1:size(old_mesh.faces, 2), [0; 0]), 1e-12, "ByRows", true, "DataScale", 1);
    if options.carry_tolerance > 0 || nargout > 3
        % the nearest translated old interface face of each new interface face
        F_new = reshape(new_mesh.intfaces, [], 1);
        F_old = reshape(old_mesh.intfaces, [], 1);
        E_new = ends(new_mesh, F_new, options.shift);
        E_old = ends(old_mesh, F_old, [0; 0]);
        gap = max(sqrt((E_new(:, 1) - E_old(:, 1).').^2 + (E_new(:, 2) - E_old(:, 2).').^2), ...
            sqrt((E_new(:, 3) - E_old(:, 3).').^2 + (E_new(:, 4) - E_old(:, 4).').^2));
        [gap, nearest] = min(gap, [], 2);
        relative = gap ./ reshape(new_mesh.area(F_new), [], 1);
        if options.carry_tolerance > 0
            loose = ~kept_f(F_new) & ~carried(F_new) & relative <= options.carry_tolerance;
            carried(F_new(loose)) = true;
            from_c(F_new(loose)) = F_old(nearest(loose));
        end
    end
end
extended = find(~kept_f & ~carried);
owner = zeros(n_faces, 1);
owner(new_mesh.boundaries(:)) = repelem(1:n_new, 3);
if ~isempty(extended)
    % k + 2 Gauss points, exact to degree 2k + 3: the exact L2 projection
    % onto Legendre P_k of a trace of degree k + 1 at most (u*)
    if options.slid_subdivisions == 1 || options.slid_subdivisions == Inf
        [xyf, tf, wf] = quadpoints_faces_2d(new_mesh, k + 1, extended, output = "cell");
    else
        [xyf, tf, wf] = composite_face_points(new_mesh, k + 1, extended, options.slid_subdivisions);
    end
    lw = wf .* legendre_basis(k, tf);
    fit_f = (lw.'*legendre_basis(k, tf)) \ lw.';                              % d1 x nqf
    rs_f = old_mesh.get_2d_local(reshape([xyf{1}(:).'; xyf{2}(:).'], 2, numel(wf), []), ...
        source(owner(extended)));
    % the rebuilt faces on the boundary of the new carve
    sliding = accumarray(new_mesh.boundaries(:), 1, [n_faces, 1]) == 1;
end
for L = 1:3
    name = skeleton(L);
    if ~has(name), continue, end
    U_hat = s.(name);
    V = zeros(size(U_hat, 1), n_faces);
    V(:, kept_f) = U_hat(:, from_f(kept_f));
    if any(carried)
        V(:, carried) = options.untranslated.(name)(:, from_c(carried));
    end
    if ~isempty(extended)
        if L == 1
            U = U_star;
            if isempty(U), U = postprocessed_displacement(sp, old_mesh, k, s_old); end
            degree = k + 1;
        else
            U = s_old.(volume(L));
            degree = k;
        end
        V(:, extended) = trace_projection(U, degree, rs_f, source(owner(extended)), fit_f);
        if (L > 1 || options.slid_displacement) && ~isempty(options.untranslated)
            on = sliding(extended);
            if options.slid_subdivisions == Inf
                V(:, extended(on)) = slid_traces_exact(options.untranslated.(name), ...
                    options.slid_correction*options.untranslated.(volume(L)), k, old_mesh, ...
                    new_mesh, extended(on), options.shift);
            else
                V(:, extended(on)) = slid_traces(options.untranslated.(name), ...
                    options.slid_correction*options.untranslated.(volume(L)), k, old_mesh, ...
                    xyf{1}(:, on), xyf{2}(:, on), options.shift, fit_f);
            end
        end
    end
    s.(name) = V;
end

if nargout > 3
    info = struct("n_carried", nnz(carried), "mismatch", NaN);
    if ~isempty(options.untranslated) && ~isempty(relative)
        info.mismatch = max(relative);
    end
end

% sizes
for name = [volume, stress, rotation]
    if has(name)
        assert(size(s.(name), 2) == n_new, "remap_solid_state:size", ...
            "%s has %d columns for %d elements.", name, size(s.(name), 2), n_new);
    end
end
for name = skeleton
    if has(name)
        assert(size(s.(name), 2) == n_faces, "remap_solid_state:size", ...
            "%s has %d columns for %d faces.", name, size(s.(name), 2), n_faces);
    end
end
end

function V = trace_projection(U, degree, rs, elements, fit)
% the Legendre dofs (component blocks) of the traces of the Dubiner
% fields U (component blocks, degree `degree`) of the given old elements
% at the reference points rs (2 x nqf x n_faces), one face per page
d = nchoosek(degree + 2, 2);
n_c = size(U, 1)/d;
d1 = size(fit, 1);
V = zeros(n_c*d1, numel(elements));
for j = 1:numel(elements)
    Q = dubiner_basis(degree, rs(:, :, j));
    for c = 0:n_c - 1
        V(c*d1 + (1:d1), j) = fit*(Q*U(c*d + (1:d), elements(j)));
    end
end
end

function V = slid_traces(U_hat, U, k, old_mesh, x, y, shift, fit)
% the Legendre dofs (component blocks) on new boundary faces, whose
% quadrature points are the columns (x, y) (nqf x n_faces), of the old
% history at z = point - shift: the old boundary trace U_hat at the
% closest point p of the old boundary to z, plus u_K(z) - u_K(p), u_K the
% Dubiner field U (degree k) of the old element owning that face. An
% interface face that is an old face translated by shift gets the old
% trace exactly; on a curve that moved, z lies within O(h^2) of the old
% polygon and the correction is of the order of that distance times the
% error of grad u_K.
F = size(old_mesh.faces, 2);
boundary = find(accumarray(old_mesh.boundaries(:), 1, [F, 1]) == 1).';
owner = zeros(F, 1);
owner(old_mesh.boundaries(:)) = repelem(1:size(old_mesh.elements, 2), 3);
a = old_mesh.coordinates(:, old_mesh.faces(1, boundary));
ab = old_mesh.coordinates(:, old_mesh.faces(2, boundary)) - a;
d1 = size(fit, 1);
d2 = nchoosek(k + 2, 2);
z = [x(:).'; y(:).'] - shift;                                       % 2 x n
tau = min(max(((z(1, :).' - a(1, :)).*ab(1, :) + (z(2, :).' - a(2, :)).*ab(2, :)) ...
    ./ sum(ab.^2, 1), 0), 1);                                      % n x n_boundary
dist = (z(1, :).' - a(1, :) - tau.*ab(1, :)).^2 + (z(2, :).' - a(2, :) - tau.*ab(2, :)).^2;
[~, i] = min(dist, [], 2);
tau = tau(sub2ind(size(tau), (1:size(z, 2)).', i));
faces = boundary(i);
p = a(:, i) + ab(:, i).*tau.';
values = zeros(size(z, 2), 2);
l = legendre_basis(k, tau);                                        % n x d1
for q = 1:size(z, 2)
    K = owner(faces(q));
    P = dubiner_basis(k, old_mesh.get_2d_local([z(:, q), p(:, q)], K));
    for c = 0:1
        values(q, c + 1) = l(q, :)*U_hat(c*d1 + (1:d1), faces(q)) ...
            + (P(1, :) - P(2, :))*U(c*d2 + (1:d2), K);
    end
end
values = reshape(values, size(x, 1), [], 2);
V = [fit*values(:, :, 1); fit*values(:, :, 2)];
end

function V = slid_traces_exact(U_hat, U, k, old_mesh, new_mesh, faces, shift)
% slid_traces projected exactly: the new face, translated back by shift,
% z(s) = z0 + s e, is split at the s where the unclamped parameter
% tau_i(s) of the closest point on an old boundary face i crosses 0 or 1;
% on each piece the closest old face (chosen at the piece's midpoint) is
% fixed, and the slid trace is a polynomial of degree k in s, integrated
% by Gauss with k + 3 points. Exact where the old polygon is convex near
% the face; at a reflex old vertex the switch between two equally distant
% faces is not a split point
F = size(old_mesh.faces, 2);
boundary = find(accumarray(old_mesh.boundaries(:), 1, [F, 1]) == 1).';
owner = zeros(F, 1);
owner(old_mesh.boundaries(:)) = repelem(1:size(old_mesh.elements, 2), 3);
a = old_mesh.coordinates(:, old_mesh.faces(1, boundary));
ab = old_mesh.coordinates(:, old_mesh.faces(2, boundary)) - a;
len2 = sum(ab.^2, 1);
d1 = k + 1;
d2 = nchoosek(k + 2, 2);
[t0, w0] = get_1d_reference_quadrature(k + 2);
V = zeros(2*d1, numel(faces));
for j = 1:numel(faces)
    x1 = new_mesh.coordinates(:, new_mesh.faces(1, faces(j)));
    e = new_mesh.coordinates(:, new_mesh.faces(2, faces(j))) - x1;
    z0 = x1 - shift;
    alpha = ((z0(1) - a(1, :)).*ab(1, :) + (z0(2) - a(2, :)).*ab(2, :))./len2;
    beta = (e(1)*ab(1, :) + e(2)*ab(2, :))./len2;
    % the old faces within reach of the segment
    near = vecnorm(z0 + e/2 - a - ab/2, 2, 1) < norm(e) + sqrt(len2);
    breaks = [-alpha(near)./beta(near), (1 - alpha(near))./beta(near)];
    breaks = unique([0, breaks(isfinite(breaks) & breaks > 0 & breaks < 1), 1]);
    breaks = breaks([true, diff(breaks) > 1e-13]);
    breaks(end) = 1;
    M = zeros(d1);
    b = zeros(d1, 2);
    for piece = 1:numel(breaks) - 1
        s0 = breaks(piece);
        h = breaks(piece + 1) - s0;
        sq = s0 + h*t0(:);                                             % nq x 1
        wq = h*w0(:);
        zm = z0 + (s0 + h/2)*e;
        tm = min(max(((zm(1) - a(1, :)).*ab(1, :) + (zm(2) - a(2, :)).*ab(2, :))./len2, 0), 1);
        [~, i] = min((zm(1) - a(1, :) - tm.*ab(1, :)).^2 + (zm(2) - a(2, :) - tm.*ab(2, :)).^2);
        tau = min(max(alpha(i) + beta(i)*sq, 0), 1);
        z = z0 + e*sq.';                                               % 2 x nq
        p = a(:, i) + ab(:, i)*tau.';
        K = owner(boundary(i));
        nq = numel(sq);
        P = dubiner_basis(k, old_mesh.get_2d_local([z, p], K));
        l = legendre_basis(k, tau);
        values = zeros(nq, 2);
        for c = 0:1
            values(:, c + 1) = l*U_hat(c*d1 + (1:d1), boundary(i)) ...
                + (P(1:nq, :) - P(nq + 1:end, :))*U(c*d2 + (1:d2), K);
        end
        Lq = legendre_basis(k, sq);
        M = M + Lq.'*(wq.*Lq);
        b = b + Lq.'*(wq.*values);
    end
    X = M \ b;
    V(:, j) = X(:);
end
end

function P = peers_points(k, r, G)
% PEERS_basis at the reference points r (2 x n), each with its own
% inverse affine map G (2 x 2 x n): n x nb x 4 (the 2 x 2 values
% column-major); the evaluation of LinearElasticityProblem's private
% shifted_projection
n = size(r, 2);
p = dubiner_basis(k, r);
d2 = size(p, 2);
P = zeros(n, 4*d2 + (k > 0)*(k + 1), 4);
for c = 1:4
    P(:, (c - 1)*d2 + (1:d2), c) = p;
end
if k == 0
    return
end
[~, d_Px, d_Py, dd_Pxx, dd_Pyy, dd_Pxy] = DubinerBasis.evaluate_basis((d2 - k):d2, r(1, :).', r(2, :).');
[b, DBx, DBy] = bubble_eval(r(1, :).', r(2, :).');
g = @(i, j) reshape(G(i, j, :), n, 1);
DP  = @(i) g(1, i).*d_Px + g(2, i).*d_Py;
DB  = @(i) g(1, i).*DBx + g(2, i).*DBy;
DDP = @(i, j) dd_Pxx.*(g(1, i).*g(1, j)) + dd_Pxy.*(g(2, i).*g(1, j) + g(1, i).*g(2, j)) ...
    + dd_Pyy.*(g(2, i).*g(2, j));
bubble = 4*d2 + (1:k + 1);
P(:, bubble, 1) =   DDP(1, 2).*b + DP(1).*DB(2);
P(:, bubble, 2) =   DDP(2, 2).*b + DP(2).*DB(2);
P(:, bubble, 3) = -(DDP(1, 1).*b + DP(1).*DB(1));
P(:, bubble, 4) = -(DDP(2, 1).*b + DP(2).*DB(1));
end

function c = centroid(m)
c = squeeze(mean(reshape(m.coordinates(:, m.elements), 2, 3, []), 2)).';
end

function fit = fit_matrix(P, w)
% values at the quadrature points -> L2-projection dofs
fit = (P.'*(w(:).*P)) \ (P.'.*w(:).');
end

function E = remap_dubiner(D, k, d2, fit, kept, source, moved, rs_old, n_new)
E = zeros(size(D, 1), n_new);
E(:, kept) = D(:, source(kept));
for j = 1:numel(moved)
    Q = dubiner_basis(k, rs_old(:, :, j));
    for c = 0:size(D, 1)/d2 - 1
        rows = c*d2 + (1:d2);
        E(rows, moved(j)) = fit*(Q*D(rows, source(moved(j))));
    end
end
end

function [xy, t, w] = composite_face_points(mesh, order, faces, m)
% the rule get_1d_reference_quadrature(order) on m equal subintervals of
% each face, in the output form of quadpoints_faces_2d (output = "cell")
[t0, w0] = get_1d_reference_quadrature(order);
t = reshape((t0(:) + (0:m - 1))/m, 1, []);
w = repmat(w0(:)/m, m, 1);
x1 = mesh.coordinates(:, mesh.faces(1, faces));
F = mesh.coordinates(:, mesh.faces(2, faces)) - x1;
xy = {x1(1, :) + t(:).*F(1, :); x1(2, :) + t(:).*F(2, :)};
end
