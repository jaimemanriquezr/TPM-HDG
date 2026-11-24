function [s, switched, moved, flux] = transport_solid_state(sp, s, old_mesh, new_mesh, k, old_keys, new_keys, ...
    old_face_keys, new_face_keys, foot, grad_foot, options)
%TRANSPORT_SOLID_STATE  Newmark solid history carried along the feet of a step, from one carve to another.
%
%   [s, switched, moved, flux] = transport_solid_state(sp, s, old_mesh, new_mesh, k, ...
%       old_keys, new_keys, old_face_keys, new_face_keys, foot, grad_foot, ...
%       quadrature_degree = 2k + 2, foot_integration = "clipped")
%
%   The semi-Lagrangian form of LinearElasticityProblem.translate_state
%   followed by remap_solid_state. s is a Newmark state.solid on old_mesh
%   (LinearElasticityProblem sp, degree k); foot maps the points x of the
%   new configuration (2 x n) to their feet y (2 x n), the positions at the
%   old time of the material points at x, and grad_foot maps them to the
%   Jacobian of the foot map (2 x 2 x n, entry (i, j) = d foot_i / d x_j).
%   Elements and faces are identified across the carves by their keys, as
%   in remap_solid_state: an element of the new carve present in the old
%   one with the same vertices is kept; the others (switched on or
%   reshaped) are moved. A foot is evaluated with the old element that
%   holds it, or, outside the old carve, with the old element with the
%   nearest centroid, extrapolated (locate_feet). Per level (displacement,
%   velocity, acceleration), on each new element K, with quadrature points
%   x_q of K:
%     velocity, acceleration (Dubiner P_k)  Pi_k[w(foot(x))], the fixed-
%       frame components unrotated;
%     displacement  on a kept element u_h|K + Pi_k[u*(foot(x)) - u*_K(x)],
%       elsewhere Pi_k[u*(foot(x))], u* the P_k+1 postprocessing
%       (postprocessed_displacement): the increment form of translate_state;
%     stress and rotation  the gradient of the level's transported volume
%       field: with G = c1 sigma - c2 tr(sigma) I + gamma the reconstructed
%       gradient (row c is grad of component c, gamma the scalar times
%       [0 1; -1 0]), G~(x) = G(foot(x)) grad_foot(x) by the chain rule;
%       the stress is the L2 projection onto PEERS(K) of C sym(G~) and the
%       rotation the P_k projection of the (1, 2) entry of skew(G~). The
%       Newmark recursions of stress and rotation then difference the
%       gradients of the transported volume levels, i.e. they return grad
%       of the velocity, which the kinematic row reads. For a uniform foot
%       (grad_foot = I) this is the projection of G(foot(x)) of the
%       translate;
%     skeleton (Legendre)  the displacement trace on a face present in both
%       carves with the same end points: the old trace plus its extension
%       from x to foot(x), int_0^1 G(x + s (foot(x) - x)) (foot(x) - x) ds,
%       with the level's own G in the old element owning the face
%       (extrapolated); on any other face, the Legendre projection of
%       u*(foot(x)). The velocity and acceleration traces on a boundary face
%       of the new carve: the old boundary trace at the closest point p of
%       the old boundary faces to z = foot(x), plus w_K(z) - w_K(p), w_K
%       the level's volume polynomial of the old element owning that face
%       (remap_solid_state's slid traces with x - shift replaced by the
%       foot); on an interior face present in both carves with the same end
%       points, the extension rule of the displacement with the level's own
%       G; on any other interior face, the Legendre projection of the
%       level's volume polynomial at the feet.
%   The integrals of the projections (volume, stress, rotation, and the
%   skeleton projections at the feet) follow foot_integration:
%   "clipped" (default; an affine foot on each element and face) integrates
%   the piecewise polynomial history at the feet exactly on the overlaps
%   of the image of each new element (face) with the old elements
%   (foot_quadrature, foot_face_quadrature), the part of the image outside
%   the old carve taking the source's polynomial, extrapolated (the old
%   element of the same vertices, else the one holding the foot of the
%   centroid; for a face, the one holding the foot of its midpoint);
%   "quadrature" uses a rule on each new element (face) with every foot
%   evaluated in the old element that holds it (locate_feet), exact on
%   polynomials for an affine foot and otherwise a quadrature of an
%   integrand that is discontinuous where the feet cross old element
%   faces. quadrature_degree is the degree of exactness of the element
%   rules (on each overlap, under "clipped") for the P_k sources
%   (velocity, acceleration, stress, rotation); the P_k+1 source u* uses
%   quadrature_degree + 2; the rules have positive weights
%   (collapsed_gauss_rule); the face rules have k + 2 Gauss points. The
%   extension of a trace on a face present in both carves, and the slid
%   boundary traces, use one old element's polynomial each and are
%   evaluated at the feet of the face's Gauss points under either option.
%   Every field's column count is checked against the new mesh;
%   displacement_previous must be empty (Newmark only). switched = [on,
%   off]: the background elements entering and leaving the solid; moved,
%   the new elements that switch on or change shape.
%   flux (computed only when requested): the transported momentum flux,
%   the L2 projection onto PEERS(K) of sigma(foot(x)) grad_foot(x)^-T, sigma
%   the old PEERS stress of the displacement level (s.stress, not
%   symmetrized), with the rule and test functions of the stress levels
%   (nb x n_new, not symmetric). For an affine foot with Jacobian R,
%   (div sigma)(foot(x)) = div(sigma(foot(x)) R^-T) (row-wise divergence),
%   so that flux in place of the stress gives the momentum identity at the
%   feet (LinearElasticityProblem.consistent_acceleration).
%
%See also: remap_solid_state, foot_projection, foot_quadrature, foot_face_quadrature,
%          locate_feet, eulerian_solid_march
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
    foot function_handle
    grad_foot function_handle
    options.quadrature_degree double {mustBeScalarOrEmpty} = []
    options.foot_integration (1,1) string {mustBeMember(options.foot_integration, ...
        ["clipped", "quadrature"])} = "clipped"
end
volume   = ["displacement", "velocity", "acceleration"];
skeleton = volume + "_skeleton";
stress   = ["stress", "stress_velocity", "stress_acceleration"];
rotation = ["rotation", "rotation_velocity", "rotation_acceleration"];
known = [volume, skeleton, stress, rotation, "displacement_previous"];
unknown = setdiff(string(fieldnames(s)).', known);
assert(isempty(unknown), "transport_solid_state:unknownField", "Unknown solid state field(s): %s.", ...
    strjoin(unknown, ", "));
assert(~isfield(s, "displacement_previous") || isempty(s.displacement_previous), ...
    "transport_solid_state:scheme", "displacement_previous is set: Newmark states only.");
for name = rotation
    if isfield(s, name) && iscell(s.(name)), s.(name) = cat(1, s.(name){:}); end
end
has = @(name) isfield(s, name) && isnumeric(s.(name)) && ~isempty(s.(name)) && ~isequal(s.(name), 0);
q = options.quadrature_degree;
if isempty(q), q = 2*k + 2; end
s_old = s;
[lambda, mu] = LinearElasticityProblem.get_lame_parameters(sp.Data.young_modulus, sp.Data.poisson_ratio);
c1 = 1/(2*mu);
c2 = lambda/(2*mu*(2*lambda + 2*mu));
d2 = nchoosek(k + 2, 2);
d2p = nchoosek(k + 3, 2);

% elements: kept (same key, same vertices) and moved
n_new = size(new_mesh.elements, 2);
[present, from] = ismembertol(new_keys, old_keys, 1e-10, "ByRows", true);
switched = [nnz(~present), size(old_mesh.elements, 2) - nnz(present)];
kept = present;
kept(present) = all(reshape(abs(new_mesh.coordinates(:, new_mesh.elements(:, present)) ...
    - old_mesh.coordinates(:, old_mesh.elements(:, from(present)))) < 1e-14, 6, []), 1).';
moved = find(~kept);

% volume levels
project = @(U, degree, elements, quad) foot_projection(old_mesh, new_mesh, U, degree, k, foot, ...
    elements = elements, quadrature_degree = quad, foot_integration = options.foot_integration);
all_elements = (1:n_new).';
for name = volume(2:3)
    if has(name), s.(name) = project(s_old.(name), k, all_elements, q); end
end
U_star = [];
if has("displacement") || has("displacement_skeleton")
    U_star = postprocessed_displacement(sp, old_mesh, k, s_old);
end
if has("displacement")
    U = project(U_star, k + 1, all_elements, q + 2);
    K = find(kept);
    % the P_k projection of u*_K on its own element: its first d2 modes
    % (Dubiner is orthogonal and hierarchical)
    own = U_star([1:d2, d2p + (1:d2)], from(K));
    U(:, K) = s_old.displacement(:, from(K)) + U(:, K) - own;
    s.displacement = U;
end

% stress and rotation levels: the pulled-back gradient of each level
rule = foot_quadrature(old_mesh, new_mesh, all_elements, foot, q, options.foot_integration);
n_points = numel(rule.w);
J = grad_foot(rule.x);
if size(J, 3) == 1, J = repmat(J, 1, 1, n_points); end
[rs, w] = collapsed_gauss_rule(q);
nq = numel(w);
G_new = get_affine_map(new_mesh, all_elements, inverse = true);
Phi = PEERS_basis(k, rs(1, :), rs(2, :), G_new);                       % nq x nb x 2 x 2 x n_new
nb = size(Phi, 2);
Phi = reshape(permute(Phi, [1 3 4 2 5]), 4*nq, nb, n_new);             % (q, i, j) x nb x K
W = repmat(w(:), 4, 1);
mass = pagemtimes(Phi, 'transpose', W.*Phi, 'none');                   % nb x nb x K
scale = 1./sqrt(reshape(mass(repmat(logical(eye(nb)), 1, 1, n_new)), nb, 1, n_new));
P = dubiner_basis(k, rs);
mass_d = P.'*(w(:).*P);
test_d = rule.w .* dubiner_basis(k, rule.rx);                          % n_points x d2
per_element = sparse(rule.element, 1:n_points, 1, n_new, n_points);
% the weighted PEERS test functions at the rule's points of each element,
% (p, i, j) x nb
at = accumarray(rule.element(:), (1:n_points).', [n_new, 1], @(v) {sort(v)});
test_s = cell(n_new, 1);
for K = 1:n_new
    idx = at{K};
    test_s{K} = reshape(permute(PEERS_basis(k, rule.rx(1, idx), rule.rx(2, idx), G_new(:, :, K)), ...
        [1 3 4 2]), [], nb) .* repmat(rule.w(idx), 4, 1);
end
% the L2 projection onto PEERS of each new element of the tensor values
% S (n x (i, j)) at the rule's points
peers_projection = @(S) reshape(scale.*pagemldivide(scale.*mass.*pagetranspose(scale), ...
    scale.*element_loads(S, test_s, at, nb, n_new)), nb, n_new);
for i = 1:3
    if ~has(stress(i)), continue, end
    G = gradient_at(old_mesh, k, s_old.(stress(i)), s_old.(rotation(i)), c1, c2, rule.L, rule.r);   % 2 x 2 x n
    G = pagemtimes(G, J);
    sym_G = (G + pagetranspose(G))/2;
    tr = sym_G(1, 1, :) + sym_G(2, 2, :);
    S = reshape(2*mu*sym_G + lambda*tr.*eye(2), 4, []).';              % C sym(G~), n x (i, j)
    s.(stress(i)) = peers_projection(S);
    s.(rotation(i)) = mass_d \ (per_element*(test_d.*reshape((G(1, 2, :) - G(2, 1, :))/2, [], 1))).';
end
if nargout > 3
    % the transported momentum flux sigma(foot(x)) grad_foot(x)^-T
    F = pagemtimes(peers_at(old_mesh, k, s_old.stress, rule.L, rule.r), 'none', pageinv(J), 'transpose');
    flux = peers_projection(reshape(F, 4, []).');
end

% skeleton levels
n_faces = size(new_mesh.faces, 2);
[present_f, from_f] = ismembertol(new_face_keys, old_face_keys, 1e-10, "ByRows", true);
kept_f = present_f;
kept_f(present_f) = all(abs(new_mesh.coordinates(:, new_mesh.faces(1, present_f)) ...
    - old_mesh.coordinates(:, old_mesh.faces(1, from_f(present_f)))) < 1e-14, 1).' ...
    & all(abs(new_mesh.coordinates(:, new_mesh.faces(2, present_f)) ...
    - old_mesh.coordinates(:, old_mesh.faces(2, from_f(present_f)))) < 1e-14, 1).';
boundary = accumarray(new_mesh.boundaries(:), 1, [n_faces, 1]) == 1;
if any(arrayfun(has, skeleton))
    % k + 2 Gauss points, exact to degree 2k + 3: the exact Legendre P_k
    % projection of a trace of degree k + 1 at most (u*)
    [xyf, tf, wf] = quadpoints_faces_2d(new_mesh, k + 1, (1:n_faces).', output = "cell");
    nqf = numel(wf);
    lw = wf .* legendre_basis(k, tf);
    fit_f = (lw.'*legendre_basis(k, tf)) \ lw.';                         % d1 x nqf
    xf = [xyf{1}(:).'; xyf{2}(:).'];                                      % 2 x (nqf*n_faces)
    yf = foot(xf);
    % the projections at the feet, on the faces not present in both carves
    rule_f = foot_face_quadrature(old_mesh, new_mesh, find(~kept_f), foot, k + 2, options.foot_integration);
    % the old element owning each old face
    F_old = size(old_mesh.faces, 2);
    owner_old = zeros(F_old, 1);
    owner_old(old_mesh.boundaries(:)) = repelem(1:size(old_mesh.elements, 2), 3);
    points = @(faces) reshape((faces(:).' - 1)*nqf + (1:nqf).', 1, []);   % columns of xf
end
d1 = k + 1;
for i = 1:3
    name = skeleton(i);
    if ~has(name), continue, end
    V = zeros(2*d1, n_faces);
    if i == 1
        extended = find(kept_f);
        projected = find(~kept_f);
        source = U_star;  degree = k + 1;
    else
        extended = find(kept_f & ~boundary);
        projected = find(~kept_f & ~boundary);
        source = s_old.(volume(i));  degree = k;
    end
    % the old trace extended from x to foot(x) along the level's gradient
    if ~isempty(extended)
        V(:, extended) = s_old.(name)(:, from_f(extended)) + extension(old_mesh, k, ...
            s_old.(stress(i)), s_old.(rotation(i)), c1, c2, owner_old(from_f(extended)), ...
            xf(:, points(extended)), yf(:, points(extended)), nqf, fit_f);
    end
    % the projection of the source's polynomial at the feet
    if ~isempty(projected)
        V(:, projected) = face_projection(rule_f, find(~kept_f), projected, source, degree, k);
    end
    % the velocity and acceleration traces slid along the feet on the boundary
    if i > 1 && any(boundary)
        on = find(boundary);
        V(:, on) = slid_traces(s_old.(name), s_old.(volume(i)), k, old_mesh, yf(:, points(on)), ...
            nqf, fit_f);
    end
    s.(name) = V;
end

% sizes
for name = [volume, stress, rotation]
    if has(name)
        assert(size(s.(name), 2) == n_new, "transport_solid_state:size", ...
            "%s has %d columns for %d elements.", name, size(s.(name), 2), n_new);
    end
end
for name = skeleton
    if has(name)
        assert(size(s.(name), 2) == n_faces, "transport_solid_state:size", ...
            "%s has %d columns for %d faces.", name, size(s.(name), 2), n_faces);
    end
end
end

function G = gradient_at(mesh, k, sigma, gamma, c1, c2, L, r)
% the reconstructed gradient c1 sigma - c2 tr(sigma) I + gamma [0 1; -1 0]
% of the elements L (1 x n) at their reference points r (2 x n): 2 x 2 x n
n = numel(L);
if iscell(gamma), gamma = cat(1, gamma{:}); end
d2 = nchoosek(k + 2, 2);
A = dubiner_basis(k, r);                                               % n x d2
g = sum(A .* gamma(1:d2, L).', 2).';
G = peers_at(mesh, k, sigma, L, r);
tr = G(1, 1, :) + G(2, 2, :);
G = c1*G - c2*tr.*eye(2);
G(1, 2, :) = G(1, 2, :) + reshape(g, 1, 1, n);
G(2, 1, :) = G(2, 1, :) - reshape(g, 1, 1, n);
end

function S = peers_at(mesh, k, sigma, L, r)
% the PEERS stress sigma of the elements L (1 x n) at their reference
% points r (2 x n): 2 x 2 x n
S = zeros(2, 2, numel(L));
Ginv = get_affine_map(mesh, (1:size(mesh.elements, 2)).', inverse = true);
[elements, ~, which] = unique(L);
for e = 1:numel(elements)
    idx = find(which == e).';
    P = PEERS_basis(k, r(1, idx), r(2, idx), Ginv(:, :, elements(e)));   % m x nb x 2 x 2
    m = numel(idx);
    T = tensorprod(P, sigma(:, elements(e)), 2, 1);                    % m x 2 x 2
    S(:, :, idx) = permute(reshape(T, m, 2, 2), [2 3 1]);
end
end

function b = element_loads(S, test_s, at, nb, n_new)
% the loads of the tensor values S (n x (i, j)) against the weighted PEERS
% test functions of each element: nb x 1 x n_new
b = zeros(nb, 1, n_new);
for K = 1:n_new
    b(:, 1, K) = test_s{K}.'*reshape(S(at{K}, :), [], 1);
end
end

function V = extension(mesh, k, sigma, gamma, c1, c2, owners, x, y, nqf, fit_f)
% the Legendre dofs (component blocks) of int_0^1 G(x + s (y - x)) (y - x) ds
% on faces (nqf points each, columns of x and y), G the reconstructed
% gradient of the face's owner element (owners, one per face), extrapolated
[s, ws] = lgwt(k + 2, 0, 1);
ns = numel(s);
n = size(x, 2);
owner = repelem(owners(:).', nqf);                                     % 1 x n
p = reshape(x, 2, 1, n) + reshape(s, 1, ns) .* reshape(y - x, 2, 1, n);   % 2 x ns x n
p = reshape(p, 2, []);
L = repelem(owner, ns);
r = reshape(mesh.get_2d_local(reshape(p, 2, 1, []), L(:)), 2, []);
G = gradient_at(mesh, k, sigma, gamma, c1, c2, L, r);                  % 2 x 2 x (ns*n)
m = reshape(repmat(reshape(y - x, 2, 1, n), 1, ns, 1), 2, 1, []);
values = reshape(pagemtimes(G, m), 2, ns, n);
values = reshape(sum(values .* reshape(ws, 1, ns), 2), 2, nqf, []);    % 2 x nqf x n_faces
V = [fit_f*reshape(values(1, :, :), nqf, []); fit_f*reshape(values(2, :, :), nqf, [])];
end

function V = slid_traces(U_hat, U, k, old_mesh, z, nqf, fit)
% remap_solid_state's slid_traces with the points z (2 x (nqf*n_faces))
% given: the old boundary trace U_hat at the closest point p of the old
% boundary to z, plus u_K(z) - u_K(p), u_K the Dubiner field U (degree k) of
% the old element owning that face; Legendre dofs on the new faces
F = size(old_mesh.faces, 2);
boundary = find(accumarray(old_mesh.boundaries(:), 1, [F, 1]) == 1).';
owner = zeros(F, 1);
owner(old_mesh.boundaries(:)) = repelem(1:size(old_mesh.elements, 2), 3);
a = old_mesh.coordinates(:, old_mesh.faces(1, boundary));
ab = old_mesh.coordinates(:, old_mesh.faces(2, boundary)) - a;
d1 = size(fit, 1);
d2 = nchoosek(k + 2, 2);
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
values = reshape(values, nqf, [], 2);
V = [fit*values(:, :, 1); fit*values(:, :, 2)];
end

function V = face_projection(rule, faces, wanted, U, degree, k)
% the Legendre P_k dofs (component blocks) on the new faces wanted of the
% Dubiner field U (degree) of the old elements at the feet of the rule's
% points (foot_face_quadrature on the faces faces)
[t, w] = lgwt(k + 2, 0, 1);
l = legendre_basis(k, t);
mass = l.'*(w(:).*l);
[in, position] = ismember(faces(rule.face), wanted);
pick = find(in);
B = dubiner_basis(degree, rule.r(:, pick));
test = rule.w(pick) .* legendre_basis(k, rule.t(pick));
S = sparse(position(pick), 1:numel(pick), 1, numel(wanted), numel(pick));
dd = nchoosek(degree + 2, 2);
d1 = k + 1;
V = zeros(2*d1, numel(wanted));
for c = 0:1
    values = sum(B .* U(c*dd + (1:dd), rule.L(pick)).', 2);
    V(c*d1 + (1:d1), :) = mass \ (S*(test.*values)).';
end
end
