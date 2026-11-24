function [V, info] = shifted_projection(mesh, k, U, shift, basis, target, k_target, extension_degree)
% Element-wise L2 projection of a translated piecewise polynomial.
%
%   V = shifted_projection(mesh, k, U, shift, basis)
%   [V, info] = shifted_projection(mesh, k, U, shift, basis, target, k_target, extension_degree)
%
% With a target mesh, the projection is onto the space of degree k_target
% of target of x -> u_h(x - shift), u_h of degree k on mesh (see
% projection_between below; extension_degree, [] when absent, is its
% q). Without one:
% basis = "dubiner" (scalar P_k, nb = d2) or "peers" (the PEERS stress
% element, 2 x 2 valued); U holds m fields of that basis stacked by rows,
% (m*nb) x N, and V the dofs of the L2 projection onto the same space of
% x -> u_h(x - shift), field by field. On element K, for every test w,
%   (V_K, w)_K = (u_K(. - shift), w)_K
%              + sum_{K'} ((u_K' - u_K)(. - shift), w)_{K ∩ (K' + shift)},
% the sum over the elements K' ≠ K that meet K - shift; the part of
% K whose preimage leaves the mesh keeps K's own polynomial, extrapolated.
% The candidates K' are the elements whose bounding box, translated,
% meets K's, searched over vertex-neighbour rings of K until a further
% ring adds none. The overlaps K ∩ (K' + shift) are clipped exactly
% (convex polygons, Sutherland-Hodgman) and integrated on a fan of
% triangles with a rule exact for the products of two basis functions,
% so the projection is exact up to round-off.
if nargin > 5
    if nargin < 8, extension_degree = []; end
    [V, info] = projection_between(mesh, target, k, k_target, U, shift, basis, extension_degree);
    return
end
N = size(mesh.elements, 2);
[B, x3] = get_affine_map(mesh);                       % x = B r + x3
G = get_affine_map(mesh, inverse = true);
det_B = reshape(B(1, 1, :).*B(2, 2, :) - B(1, 2, :).*B(2, 1, :), 1, []);

% a rule exact for the products of two basis functions (degree 2k, or
% 2k + 2 with the PEERS bubbles); sum(w) = 1
switch basis
    case "dubiner"
        evaluate = @(r, L) dubiner_basis(k, r);                  % n x nb (x 1)
        [rs, w] = get_2d_reference_quadrature(k);
    case "peers"
        evaluate = @(r, L) peers_points(k, r, G(:, :, L));       % n x nb x 4
        [rs, w] = get_2d_reference_quadrature(k + 1);
end
nq = numel(w);
nb = size(evaluate(rs(:, 1), 1), 2);
m = size(U, 1)/nb;
U = reshape(U, nb, m, N);
% basis values at the points y (2 x nq x n_pages) of elements L
% (1 x n_pages), as nb x (nv*nq) x n_pages
at = @(y, L) values(evaluate, G, x3, y, L, nq);
% the moments (w_i, u_L(y))_q of the fields of elements L at the points
% y against the test functions Phi (nb x nv*nq x n_pages), with the
% weights W (1 x nq x n_pages); nb x m x n_pages
moments_of = @(Phi, W, L, y) pagemtimes(Phi .* repmat(W, 1, size(Phi, 2)/nq), 'none', ...
    pagemtimes(U(:, :, L), 'transpose', at(y, L), 'none'), 'transpose');

% the element's own polynomial, translated: reference points on K
x_K = pagemtimes(B, rs) + reshape(x3, 2, 1, N);       % 2 x nq x N
Phi_K = at(x_K, 1:N);
W_K = repmat(reshape(w, 1, nq), 1, 1, N);
mass = pagemtimes(Phi_K .* repmat(W_K, 1, size(Phi_K, 2)/nq), 'none', Phi_K, 'transpose');
moments = moments_of(Phi_K, W_K, 1:N, x_K - shift);

% the overlaps with the other elements
vertices = @(L) reshape(mesh.coordinates(1:2, mesh.elements(:, L)), 2, 3, []);
[K, Kp] = candidate_pairs(mesh, vertices(1:N), shift);
[P, n] = clip_polygons(counterclockwise(vertices(K)), counterclockwise(vertices(Kp) + shift));

% fan triangulation of each overlap: (P_1, P_j, P_j+1), j = 2..n-1
pair = [];  a = [];  b = [];  c = [];
for j = 2:size(P, 2) - 1
    has = find(n > j);
    pair = [pair, has]; %#ok<AGROW>
    a = [a, reshape(P(:, 1, has), 2, [])]; %#ok<AGROW>
    b = [b, reshape(P(:, j, has), 2, [])]; %#ok<AGROW>
    c = [c, reshape(P(:, j + 1, has), 2, [])]; %#ok<AGROW>
end
area_T = ((b(1, :) - a(1, :)).*(c(2, :) - a(2, :)) ...
        - (b(2, :) - a(2, :)).*(c(1, :) - a(1, :)))/2;
% overlaps of zero area (a shared edge or vertex) contribute nothing
keep = area_T > 1e-14*abs(det_B(K(pair)));
pair = pair(keep);  a = a(:, keep);  b = b(:, keep);  c = c(:, keep);  area_T = area_T(keep);
n_tri = numel(pair);
if n_tri > 0
    owner = K(pair);
    % weights relative to |K|, as on K itself
    W = reshape(w, 1, nq) .* reshape(area_T ./ (abs(det_B(owner))/2), 1, 1, []);
    x = reshape(a, 2, 1, []) + reshape(b - a, 2, 1, []).*rs(1, :) ...
        + reshape(c - a, 2, 1, []).*rs(2, :);                       % 2 x nq x n_tri
    Phi = at(x, owner);
    jump = moments_of(Phi, W, Kp(pair), x - shift) - moments_of(Phi, W, owner, x - shift);
    moments = moments + reshape(reshape(jump, nb*m, n_tri) ...
        * sparse(1:n_tri, owner, 1, n_tri, N), nb, m, N);
end

% solved with the mass matrices scaled to a unit diagonal: the PEERS bubble
% columns scale with G^2 and the Dubiner ones do not
D = 1./sqrt(reshape(mass(repmat(logical(eye(nb)), 1, 1, N)), nb, 1, N));
V = reshape(D.*pagemldivide(D.*mass.*pagetranspose(D), D.*moments), nb*m, N);
end

function [K, Kp] = candidate_pairs(mesh, T, shift)
% the pairs K ≠ K' whose bounding boxes meet once K' is translated by
% shift, over growing vertex-neighbour rings of K
N = size(T, 3);
lo = reshape(min(T, [], 2), 2, N);
hi = reshape(max(T, [], 2), 2, N);
slack = 1e-12*max(hi - lo, [], "all");
neighbours = sparse(mesh.elements(:), repelem(1:N, 3), 1);
neighbours = double(neighbours.'*neighbours > 0);
ring = neighbours;
n_pairs = -1;
while true
    [Kp, K] = find(ring);
    meets = all(lo(:, K) <= hi(:, Kp) + shift + slack & lo(:, Kp) + shift <= hi(:, K) + slack, 1) ...
        & (K ~= Kp).';
    if nnz(meets) == n_pairs, break; end
    n_pairs = nnz(meets);
    ring = double(ring*neighbours > 0);
end
K = K(meets).';
Kp = Kp(meets).';
end

function Phi = values(evaluate, G, x3, y, L, nq)
% basis of elements L at the physical points y, nb x (nv*nq) x n_pages
n_pages = size(y, 3);
r = pagemtimes(G(:, :, L), y - reshape(x3(:, L), 2, 1, []));
Phi = evaluate(reshape(r, 2, []), repelem(L, nq));               % n x nb x nv
nb = size(Phi, 2);
Phi = reshape(Phi, nq, n_pages, nb, []);                         % q, page, i, c
Phi = reshape(permute(Phi, [3 1 4 2]), nb, [], n_pages);         % i, (q, c), page
end

function P = peers_points(k, r, G)
% PEERS_basis at the reference points r (2 x n), each with its own
% inverse affine map G (2 x 2 x n): n x nb x 4 (the 2 x 2 values
% column-major), the scalar Dubiner modes times E11, E21, E12, E22,
% then the curl-curl bubble modes of curl_curl_A_dubiner_bubble_eval
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
P(:, bubble, 1) =   DDP(1, 2).*b + DP(1).*DB(2);     % (1, 1)
P(:, bubble, 2) =   DDP(2, 2).*b + DP(2).*DB(2);     % (2, 1)
P(:, bubble, 3) = -(DDP(1, 1).*b + DP(1).*DB(1));    % (1, 2)
P(:, bubble, 4) = -(DDP(2, 1).*b + DP(2).*DB(1));    % (2, 2)
end

function T = counterclockwise(T)
% triangles 2 x 3 x n with their vertices in counterclockwise order
signed = (T(1, 2, :) - T(1, 1, :)).*(T(2, 3, :) - T(2, 1, :)) ...
       - (T(2, 2, :) - T(2, 1, :)).*(T(1, 3, :) - T(1, 1, :));
flip = signed < 0;
T(:, [2 3], flip) = T(:, [3 2], flip);
end

function [P, n] = clip_polygons(subject, clip)
% Sutherland-Hodgman, vectorized over pages: the convex polygons
% subject ∩ clip, both counterclockwise triangles (2 x 3 x n). P is
% 2 x m x n, m >= 6 (the intersection of two triangles has at most six
% vertices; with round-off, near-collinear edges can give more: every
% pass may at most double the count, and the pages grow as needed), n the
% vertex count of each page (below 3 when they do not overlap).
n_pages = size(subject, 3);
MAX = 6;
P = zeros(2, MAX, n_pages);
P(:, 1:3, :) = subject;
n = 3*ones(1, n_pages);
pages = 1:n_pages;
for e = 1:3
    % room for every vertex of this pass: one crossing and one kept
    % vertex per input vertex at most
    MAX_in = size(P, 2);
    MAX = max(MAX_in, 2*max([n, 0]));
    a = reshape(clip(:, e, :), 2, []);
    b = reshape(clip(:, mod(e, 3) + 1, :), 2, []);
    side = @(p) (b(1, :) - a(1, :)).*(p(2, :) - a(2, :)) ...
              - (b(2, :) - a(2, :)).*(p(1, :) - a(1, :));
    Q = zeros(2, MAX, n_pages);
    n_new = zeros(1, n_pages);
    for i = 1:max([n, 0])
        active = i <= n;
        % the previous vertex, cyclically (vertex n before vertex 1)
        j_prev = (i - 1)*ones(1, n_pages);
        if i == 1, j_prev = max(n, 1); end
        current = reshape(P(:, i, :), 2, []);
        prior = P(:, (pages - 1)*MAX_in + j_prev);
        d_cur = side(current);
        d_pri = side(prior);
        in_cur = d_cur >= 0;
        in_pri = d_pri >= 0;
        crossing = active & (in_cur ~= in_pri);
        lambda = d_pri ./ (d_pri - d_cur + ~crossing);
        X = prior + (current - prior).*lambda;
        % the crossing point first, then the current vertex when inside
        n_new = n_new + crossing;
        Q = put(Q, X, n_new, crossing);
        keep = active & in_cur;
        n_new = n_new + keep;
        Q = put(Q, current, n_new, keep);
    end
    P = Q;
    n = n_new;
end
end

function Q = put(Q, X, slot, mask)
% Q(:, slot(p), p) = X(:, p) on the pages where mask holds
[~, MAX, n_pages] = size(Q);
p = find(mask);
Q = reshape(Q, 2, MAX*n_pages);
Q(:, (p - 1)*MAX + slot(p)) = X(:, p);
Q = reshape(Q, 2, MAX, n_pages);
end

function [V, info] = projection_between(old_mesh, new_mesh, k_from, k_to, U, shift, basis, q)
% The L2 projection onto the space of every element K of new_mesh of
% x -> u(x - shift), u the piecewise field U on old_mesh: Dubiner of degree
% k_from onto Dubiner of degree k_to, or PEERS of degree k onto PEERS of K
% (k_from = k_to = k; the bubbles of K through K's own affine map). On K,
% for every test w,
%   (V_K, w)_K = sum_{K'} (u_K'(. - shift), w)_{K ∩ (K' + shift)}
%              + (u_E(. - shift), w)_{K \ (old_mesh + shift)},
% the sum over the old elements K' that meet K - shift (old elements
% whose translated bounding box meets K's, found on a uniform grid of bins
% over the bounding box of new_mesh, bin side the largest bounding-box
% side of both meshes). The remainder R_K = K \ (old_mesh + shift), where
% x - shift lies in no old element, takes the polynomial of the old
% element E = E(K), extrapolated: E is the old element whose translate
% overlaps K with the largest area (ties: the lowest index); if no
% translate meets K, the old element whose translated centroid is nearest
% to K's centroid (ties: the lowest index). The remainder is integrated
% exactly as (u_E(. - shift), w)_K - sum_{K'} (u_E(. - shift), w)_{K ∩
% (K' + shift)}, on the elements whose relative uncovered area exceeds
% 1e-12 (elsewhere the two terms cancel to round-off). The old elements
% are assumed not to overlap one another. Overlaps are clipped and
% integrated as in the one-mesh case, so the projection is exact up to
% round-off. info.source is E (1 x N_new), info.uncovered |R_K|/|K|.
% With q in {0, 1} (Dubiner only; [] or q >= k_from: the polynomial of E),
% R_K takes instead the Taylor polynomial of degree q of u_E about x_K,
% the point of the closed triangle E + shift nearest to the centroid of
% R_K, (|K| c_K - sum_K' |K ∩ (K' + shift)| c_(K ∩ (K' + shift)))/|R_K|;
% in E's Dubiner dofs it is T_K u_E (taylor_matrices), integrated as u_E
% above. info.taylor_point is x_K (2 x N_new, NaN off the remainder).
N = size(new_mesh.elements, 2);
N_old = size(old_mesh.elements, 2);
[B, x3] = get_affine_map(new_mesh);
G = get_affine_map(new_mesh, inverse = true);
[~, x3_old] = get_affine_map(old_mesh);
G_old = get_affine_map(old_mesh, inverse = true);
det_B = reshape(B(1, 1, :).*B(2, 2, :) - B(1, 2, :).*B(2, 1, :), 1, []);

% a rule exact for the products of a test and a source function, and of
% two test functions (degree <= 2*max(k_from, k_to), or 2k + 2 with the
% PEERS bubbles); sum(w) = 1
switch basis
    case "dubiner"
        test = @(r, L) dubiner_basis(k_to, r);
        source = @(r, L) dubiner_basis(k_from, r);
        [rs, w] = get_2d_reference_quadrature(max(k_from, k_to));
    case "peers"
        assert(k_from == k_to, "LinearElasticityProblem:peersDegree", ...
            "The PEERS projection between meshes keeps the degree.");
        test = @(r, L) peers_points(k_to, r, G(:, :, L));
        source = @(r, L) peers_points(k_from, r, G_old(:, :, L));
        [rs, w] = get_2d_reference_quadrature(k_to + 1);
end
nq = numel(w);
nb = size(test(rs(:, 1), 1), 2);
nb_from = size(source(rs(:, 1), 1), 2);
m = size(U, 1)/nb_from;
U = reshape(U, nb_from, m, N_old);
% test functions of new elements L, source functions of old elements L, at
% the points y (2 x nq x n_pages), as nb x (nv*nq) x n_pages
at = @(y, L) values(test, G, x3, y, L, nq);
from = @(y, L) values(source, G_old, x3_old, y, L, nq);
% the moments (w_i, u_L(y))_q of the old fields of elements L at the points
% y against the test functions Phi, with the weights W (1 x nq x n_pages)
moments_of = @(Phi, W, L, y) pagemtimes(Phi .* repmat(W, 1, size(Phi, 2)/nq), 'none', ...
    pagemtimes(U(:, :, L), 'transpose', from(y, L), 'none'), 'transpose');

x_K = pagemtimes(B, rs) + reshape(x3, 2, 1, N);       % 2 x nq x N
Phi_K = at(x_K, 1:N);
W_K = repmat(reshape(w, 1, nq), 1, 1, N);
mass = pagemtimes(Phi_K .* repmat(W_K, 1, size(Phi_K, 2)/nq), 'none', Phi_K, 'transpose');
moments = zeros(nb, m, N);

% the overlaps K ∩ (K' + shift)
T = reshape(new_mesh.coordinates(1:2, new_mesh.elements), 2, 3, N);
T_old = reshape(old_mesh.coordinates(1:2, old_mesh.elements), 2, 3, N_old) + shift;
[K, Kp] = overlapping_boxes(T, T_old);
[P, n] = clip_polygons(counterclockwise(T(:, :, K)), counterclockwise(T_old(:, :, Kp)));
pair = zeros(1, 0);  a = zeros(2, 0);  b = a;  c = a;
for j = 2:size(P, 2) - 1
    has = find(n > j);
    pair = [pair, has]; %#ok<AGROW>
    a = [a, reshape(P(:, 1, has), 2, [])]; %#ok<AGROW>
    b = [b, reshape(P(:, j, has), 2, [])]; %#ok<AGROW>
    c = [c, reshape(P(:, j + 1, has), 2, [])]; %#ok<AGROW>
end
area_T = ((b(1, :) - a(1, :)).*(c(2, :) - a(2, :)) ...
        - (b(2, :) - a(2, :)).*(c(1, :) - a(1, :)))/2;
keep = area_T > 1e-14*abs(det_B(K(pair)));
pair = pair(keep);  a = a(:, keep);  b = b(:, keep);  c = c(:, keep);  area_T = area_T(keep);
n_tri = numel(pair);
overlap = accumarray(pair(:), area_T(:), [numel(K), 1]).';

% the uncovered fraction and the extrapolation source E
area_K = abs(det_B)/2;
uncovered = max(area_K - accumarray(K(:), overlap(:), [N, 1]).', 0)./area_K;
E = zeros(1, N);
met = overlap > 0;
if any(met)
    [~, order] = sortrows([K(met).', -overlap(met).', Kp(met).']);
    Km = K(met);  Kpm = Kp(met);
    [owners, first] = unique(Km(order), "first");
    E(owners) = Kpm(order(first));
end
alone = find(E == 0);
if ~isempty(alone)
    E(alone) = nearest_centroid(T(:, :, alone), T_old);
end
R = uncovered > 1e-12;

% the moments of the remainder's field of new elements L at the points y
remainder_of = @(Phi, W, L, y) moments_of(Phi, W, E(L), y);
if ~isempty(q)
    assert(basis == "dubiner" && any(q == [0 1]), "LinearElasticityProblem:extensionDegree", ...
        "The Taylor extension takes the Dubiner basis and the degree 0 or 1.");
end
taylor = ~isempty(q) && q < k_from;
if taylor
    x_T = nan(2, N);
    L_R = find(R);
    if ~isempty(L_R)
        % the centroid of R_K, from K's and the overlaps' first moments
        covered = zeros(2, N);
        if n_tri > 0
            centroid_T = (a + b + c)/3;
            covered = [accumarray(K(pair).', (area_T.*centroid_T(1, :)).', [N, 1]).'; ...
                       accumarray(K(pair).', (area_T.*centroid_T(2, :)).', [N, 1]).'];
        end
        covered_area = accumarray(K(:), overlap(:), [N, 1]).';
        c_R = (area_K.*reshape(mean(T, 2), 2, N) - covered)./(area_K - covered_area);
        x_T(:, L_R) = nearest_point(T_old(:, :, E(L_R)), c_R(:, L_R));
        r0 = pagemtimes(G_old(:, :, E(L_R)), reshape(x_T(:, L_R) - shift - x3_old(:, E(L_R)), 2, 1, []));
        C = zeros(nb_from, m, N);
        C(:, :, L_R) = pagemtimes(taylor_matrices(k_from, q, reshape(r0, 2, [])), U(:, :, E(L_R)));
        remainder_of = @(Phi, W, L, y) pagemtimes(Phi .* repmat(W, 1, size(Phi, 2)/nq), 'none', ...
            pagemtimes(C(:, :, L), 'transpose', from(y, E(L)), 'none'), 'transpose');
    end
end

if any(R)
    moments(:, :, R) = remainder_of(Phi_K(:, :, R), W_K(:, :, R), find(R), x_K(:, :, R) - shift);
end
if n_tri > 0
    owner = K(pair);
    W = reshape(w, 1, nq) .* reshape(area_T ./ area_K(owner), 1, 1, []);
    x = reshape(a, 2, 1, []) + reshape(b - a, 2, 1, []).*rs(1, :) ...
        + reshape(c - a, 2, 1, []).*rs(2, :);                       % 2 x nq x n_tri
    Phi = at(x, owner);
    piece = moments_of(Phi, W, Kp(pair), x - shift);
    in_R = find(R(owner));
    if ~isempty(in_R)
        piece(:, :, in_R) = piece(:, :, in_R) ...
            - remainder_of(Phi(:, :, in_R), W(:, :, in_R), owner(in_R), x(:, :, in_R) - shift);
    end
    moments = moments + reshape(reshape(piece, nb*m, n_tri) ...
        * sparse(1:n_tri, owner, 1, n_tri, N), nb, m, N);
end

D = 1./sqrt(reshape(mass(repmat(logical(eye(nb)), 1, 1, N)), nb, 1, N));
V = reshape(D.*pagemldivide(D.*mass.*pagetranspose(D), D.*moments), nb*m, N);
info = struct("source", E, "uncovered", uncovered);
if taylor
    info.taylor_point = x_T;
end
end

function T = taylor_matrices(k, q, r0)
% the Dubiner dofs (degree k) of the Taylor polynomials of degree q in
% {0, 1} of the Dubiner basis functions about the reference points r0
% (2 x n): T(:, i, p) expands the i-th function about r0(:, p), by an L2
% fit on the reference triangle (exact, q <= k); nb x nb x n
[rs, w] = get_2d_reference_quadrature(k);
Phi = dubiner_basis(k, rs);                                     % nq x nb
[P0, dPr, dPs] = dubiner_basis(k, r0);                          % n x nb
[nb, n] = deal(size(Phi, 2), size(r0, 2));
values = repmat(reshape(P0.', 1, nb, n), numel(w), 1, 1);       % nq x nb x n
if q == 1
    values = values + reshape(dPr.', 1, nb, n).*(rs(1, :).' - reshape(r0(1, :), 1, 1, n)) ...
        + reshape(dPs.', 1, nb, n).*(rs(2, :).' - reshape(r0(2, :), 1, 1, n));
end
fit = (Phi.'*(w(:).*Phi)) \ (Phi.'.*w(:).');                   % nb x nq
T = pagemtimes(fit, values);
end

function p = nearest_point(T, x)
% the point of the closed triangle T(:, :, i) nearest to x(:, i), 2 x n
n = size(x, 2);
T = counterclockwise(T);
p = x;
best = inf(1, n);
inside = true(1, n);
for e = 1:3
    a = reshape(T(:, e, :), 2, n);
    d = reshape(T(:, mod(e, 3) + 1, :), 2, n) - a;
    inside = inside & d(1, :).*(x(2, :) - a(2, :)) - d(2, :).*(x(1, :) - a(1, :)) >= 0;
    y = a + min(max(sum((x - a).*d, 1)./sum(d.^2, 1), 0), 1).*d;
    dist = sum((x - y).^2, 1);
    closer = dist < best;
    best(closer) = dist(closer);
    p(:, closer) = y(:, closer);
end
p(:, inside) = x(:, inside);
end

function [K, Kp] = overlapping_boxes(T, S)
% the pairs (K, K') of triangles T(:, :, K), S(:, :, K') whose bounding
% boxes meet, through a uniform grid of bins over T's bounding box with the
% largest box side of T and S as the bin side
lo = reshape(min(T, [], 2), 2, []);   hi = reshape(max(T, [], 2), 2, []);
lo_s = reshape(min(S, [], 2), 2, []); hi_s = reshape(max(S, [], 2), 2, []);
width = max([hi - lo, hi_s - lo_s], [], "all");
slack = 1e-12*width;
origin = min(lo, [], 2) - slack;
top = max(hi, [], 2) + slack;
n_bins = max(ceil((top - origin)/width), 1);
near = find(all(hi_s + slack >= origin & lo_s - slack <= top, 1));
A = bin_incidence(lo - slack, hi + slack, origin, width, n_bins);
A_s = bin_incidence(lo_s(:, near) - slack, hi_s(:, near) + slack, origin, width, n_bins);
[K, j] = find(A*A_s.');
Kp = near(j);
K = reshape(K, 1, []);
Kp = reshape(Kp, 1, []);
meets = all(lo(:, K) <= hi_s(:, Kp) + slack & lo_s(:, Kp) <= hi(:, K) + slack, 1);
K = K(meets);
Kp = Kp(meets);
end

function A = bin_incidence(lo, hi, origin, width, n_bins)
% boxes [lo, hi] (2 x n) x bins, the bins they meet (clamped to the grid)
bin = @(x) min(max(floor((x - origin)/width) + 1, 1), n_bins);
i0 = bin(lo);
i1 = bin(hi);
rows = [];  cols = [];
for dx = 0:max(i1(1, :) - i0(1, :))
    for dy = 0:max(i1(2, :) - i0(2, :))
        has = find(i0(1, :) + dx <= i1(1, :) & i0(2, :) + dy <= i1(2, :));
        rows = [rows, has]; %#ok<AGROW>
        cols = [cols, i0(1, has) + dx + (i0(2, has) + dy - 1)*n_bins(1)]; %#ok<AGROW>
    end
end
A = sparse(rows, cols, 1, size(lo, 2), prod(n_bins));
end

function E = nearest_centroid(T, S)
% for each triangle T(:, :, i), the triangle of S with the nearest
% centroid (ties: the lowest index), in blocks of 256
c = reshape(mean(T, 2), 2, []);
c_s = reshape(mean(S, 2), 2, []);
E = zeros(1, size(c, 2));
for i0 = 1:256:size(c, 2)
    i = i0:min(i0 + 255, size(c, 2));
    [~, E(i)] = min((c_s(1, :).' - c(1, i)).^2 + (c_s(2, :).' - c(2, i)).^2, [], 1);
end
end
