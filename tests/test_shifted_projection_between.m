% test_shifted_projection_between.m
% LinearElasticityProblem.translated_field_between: the L2 projection onto
% the elements of a new mesh of x -> u(x - s), u piecewise polynomial on an
% old mesh (Dubiner P_kf -> P_kt, or PEERS -> PEERS of the new element),
% overlaps clipped exactly, the uncovered part of a new element K (x - s in
% no old element) given the extrapolated polynomial of the old element
% info.source(K). Carves: StructuredSimplicialGrid(N, N, "crossed").split
% (circle of radius r, max_normal_deviation = 90), snapped with
% snap_to_curve(..., normal_sign = -1) unless said otherwise.
%   a. old = new mesh (carve N = 16, centre (0.5, 0.5)): equals the one-mesh
%      projection (translated_field for Dubiner; the stress of
%      translate_state for PEERS) to 1e-12 relative, k = 1, 2, for shifts
%      of 0.05 h to 0.3 h in several directions, on every element whose
%      uncovered part is empty or whose source is the element itself (the
%      one-mesh projection extrapolates the element's own polynomial); the
%      number of other elements is printed (0 expected for these shifts).
%   b. zero shift, carves at centres c1 = (0.5, 0.5) and c2 = (0.537,
%      0.511), N = 16: a random global polynomial of degree k is reproduced
%      on the new carve (Dubiner), a random P_k tensor (PEERS), a degree
%      k + 1 polynomial from P_k+1 onto P_k+1, and onto P_k it equals the
%      element-wise L2 projection on the new carve (1e-11 relative).
%   c. shift s = c2 - c1, same carves: polynomial reproduction as in b
%      (u(x - s)); for a random piecewise field the result agrees with a
%      composite-quadrature projection (nsub^2 uniform subtriangles of
%      each new element, the same extrapolation source) to that quadrature's
%      accuracy: the differences of the dofs (Frobenius norm relative to
%      the result's; also the largest, relative to the largest dof) for
%      nsub = 4, 8, 16, 32 are printed, and the first must decrease by a
%      factor >= 3 from nsub = 4 to 32 (~1/nsub: the subtriangles cut by
%      the old element edges).
%   d. conservation: the integral of the projection over the new carve
%      equals the integral of the shifted field: in case c against the
%      composite quadrature (nsub = 16, 32, 64, printed, decreasing); and to
%      1e-12 relative when the new mesh lies in old + s: old = background
%      grid N = 8 (unit square), new = the unsnapped carve of the grid
%      N = 16 (nested in it) translated cell-wise, s = (1/8, -2/8), where the
%      integral is computed per new element from the old element containing
%      it, without clipping.
%   f. exact reference: new mesh = the triangle (0,0), (1,0), (0,1), old =
%      its uniform refinement into 4^2 triangles, s = (1/16, 0): every
%      subtriangle of the composite rule with nsub = 16 lies in one old
%      element or in the uncovered strip, so that rule is exact; Dubiner
%      and PEERS, k = 1, 2, agree to 1e-10 relative.
%   g. round-off in the clipping: two triangles whose clip in floating
%      point has seven vertices (near-collinear edges, perturbations of
%      1e-16): the uncovered fraction equals the exact 3/4 to 1e-12 and
%      the result is finite.
%   e. L2 contraction ||P u(. - s)||_new <= ||u||_old (1 + 1e-12) when the
%      new mesh lies in old + s: the nested case of d, and a snapped carve
%      of radius 0.2 at c2 from a snapped carve of radius 0.25 at c1 (its
%      uncovered fractions checked to vanish).
clc; clear;
n_fail = 0;
tol = 1e-12;
c1 = [0.5; 0.5];
c2 = [0.537; 0.511];
sp = LinearElasticityProblem(2, time_scheme = "newmark");
rng(11);

%% a. old = new mesh
mesh = carve(16, c1, 0.2, true);
h = 1/16;
directions = [1 0; 0.6 0.8; -0.8 0.6; -0.28 -0.96].';
for k = 1:2
    d2 = nchoosek(k + 2, 2);
    for scale = [0.05, 0.3]
        for j = 1:size(directions, 2)
            s = scale*h*directions(:, j);
            U = randn(2*d2, size(mesh.elements, 2));
            V1 = translated_field(sp, mesh, k, U, s);
            [V2, info] = translated_field_between(sp, mesh, mesh, k, k, U, s);
            same = info.uncovered <= 1e-12 | info.source == 1:size(mesh.elements, 2);
            d = max(abs(V1(:, same) - V2(:, same)), [], "all")/max(abs(V1), [], "all");
            % PEERS: the stress level of translate_state
            state = random_state(mesh, k);
            S = translate_state(sp, mesh, state, s, k);
            [W2, info_p] = translated_field_between(sp, mesh, mesh, k, k, state.solid.stress, s, ...
                basis = "peers");
            same_p = info_p.uncovered <= 1e-12 | info_p.source == 1:size(mesh.elements, 2);
            dp = max(abs(S.solid.stress(:, same_p) - W2(:, same_p)), [], "all") ...
                /max(abs(S.solid.stress), [], "all");
            fprintf("a. k = %d s = %.3f h (%5.2f, %5.2f): dubiner %.1e, peers %.1e, other elements %d / %d\n", ...
                k, scale, directions(:, j), d, dp, nnz(~same), nnz(~same_p));
            n_fail = n_fail + check(d <= tol && dp <= tol, sprintf("a. k = %d shift %d", k, j));
        end
    end
end

%% b, c. two carves, zero shift and shift c2 - c1
old = carve(16, c1, 0.2, true);
new = carve(16, c2, 0.2, true);
for case_ = ["b", "c"]
    s = [0; 0];
    if case_ == "c", s = c2 - c1; end
    for k = 1:2
        p = random_polynomial(k, 2);
        V = translated_field_between(sp, old, new, k, k, dubiner_dofs(old, k, p, [0; 0]), s);
        d = rel(V, dubiner_dofs(new, k, p, s));
        q = random_polynomial(k, 4);
        [W, info] = translated_field_between(sp, old, new, k, k, peers_dofs(old, k, q, [0; 0]), s, ...
            basis = "peers");
        dp = rel(W, peers_dofs(new, k, q, s));
        p1 = random_polynomial(k + 1, 2);
        U1 = dubiner_dofs(old, k + 1, p1, [0; 0]);
        d11 = rel(translated_field_between(sp, old, new, k + 1, k + 1, U1, s), dubiner_dofs(new, k + 1, p1, s));
        d10 = rel(translated_field_between(sp, old, new, k + 1, k, U1, s), dubiner_dofs(new, k, p1, s));
        fprintf("%s. k = %d: P_k %.1e, PEERS P_k %.1e, P_k+1 -> P_k+1 %.1e, P_k+1 -> P_k %.1e; " + ...
            "uncovered elements %d (max fraction %.1e)\n", case_, k, d, dp, d11, d10, ...
            nnz(info.uncovered > 1e-12), max(info.uncovered));
        n_fail = n_fail + check(max([d, dp, d11, d10]) <= 1e-11, sprintf("%s. reproduction k = %d", case_, k));
    end
end

%% c, d. piecewise field against composite quadrature
s = c2 - c1;
for k = 1:2
    for basis = ["dubiner", "peers"]
        if basis == "dubiner"
            U = randn(2*nchoosek(k + 2, 2), size(old.elements, 2));
        else
            U = randn(4*nchoosek(k + 2, 2) + k + 1, size(old.elements, 2));
        end
        [V, info] = translated_field_between(sp, old, new, k, k, U, s, basis = basis);
        I = integral_of(new, k, V, basis);
        nsubs = [4, 8, 16, 32];
        if basis == "dubiner" && k == 1, nsubs = [nsubs, 64]; end %#ok<AGROW>
        [dV, dV_max, dI] = deal(zeros(size(nsubs)));
        for i = 1:numel(nsubs)
            [Vc, Ic] = composite_projection(old, new, k, k, U, s, basis, nsubs(i), info.source);
            dV(i) = norm(V - Vc, "fro")/norm(V, "fro");
            dV_max(i) = rel(V, Vc);
            dI(i) = max(abs(I - Ic))/max(abs(I));
        end
        fprintf("c. k = %d %-7s nsub = %s: dofs rms %s, max %s\n", k, basis, mat2str(nsubs), ...
            num2str(dV, " %.2e"), num2str(dV_max, " %.2e"));
        fprintf("d. k = %d %-7s nsub = %s: integral %s\n", k, basis, mat2str(nsubs), num2str(dI, " %.2e"));
        n_fail = n_fail + check(dV(4) <= dV(1)/3 && dI(end) < dI(1), ...
            sprintf("c/d. composite convergence k = %d %s", k, basis));
    end
end

%% f. exact reference: composite rule aligned with the old elements
m = 4;
[i, j] = meshgrid(0:m);
ij = [i(i + j <= m), j(i + j <= m)].';
X = ij/m;
id = @(a, b) find(ij(1, :) == a & ij(2, :) == b);
T = zeros(3, 0);
for a = 0:m - 1
    for b = 0:m - 1 - a
        T(:, end + 1) = [id(a, b); id(a + 1, b); id(a, b + 1)]; %#ok<SAGROW>
        if a + b <= m - 2
            T(:, end + 1) = [id(a + 1, b); id(a + 1, b + 1); id(a, b + 1)]; %#ok<SAGROW>
        end
    end
end
fine = SimplicialGrid(X, T, mesh_generator = "array");
single = SimplicialGrid([0 1 0; 0 0 1], [1; 2; 3], mesh_generator = "array");
s = [1/16; 0];
for k = 1:2
    for basis = ["dubiner", "peers"]
        nb = 2*nchoosek(k + 2, 2);
        if basis == "peers", nb = 4*nchoosek(k + 2, 2) + k + 1; end
        U = randn(nb, size(fine.elements, 2));
        [V, info] = translated_field_between(sp, fine, single, k, k, U, s, basis = basis);
        Vc = composite_projection(fine, single, k, k, U, s, basis, 16, info.source);
        fprintf("f. aligned  k = %d %-7s: dofs %.1e (uncovered %.4f)\n", k, basis, rel(V, Vc), ...
            info.uncovered);
        n_fail = n_fail + check(rel(V, Vc) <= 1e-10, sprintf("f. aligned k = %d %s", k, basis));
    end
end

%% g. round-off in the clipping: a seven-vertex clip
near = SimplicialGrid([1 0.5 0.25; 0.75 0.75 0], [1; 2; 3], mesh_generator = "array");
skew = SimplicialGrid([0.75, -1.2150537312839183e-16, 0.25; ...
    0.7500000000000001, -2.8192914384307714e-16, 2.6196309010184792e-17], [1; 2; 3], ...
    mesh_generator = "array");
[V, info] = translated_field_between(sp, skew, near, 1, 1, randn(6, 1), [0; 0]);
fprintf("g. seven-vertex clip: uncovered %.15f (exact 0.75)\n", info.uncovered);
n_fail = n_fail + check(abs(info.uncovered - 0.75) <= 1e-12 && all(isfinite(V(:))), ...
    "g. seven-vertex clip");

%% d, e. new mesh inside old + s
coarse = StructuredSimplicialGrid(8, 8, "crossed");
nested = carve(16, c1, 0.2, false);
s = [1/8; -2/8];
nested = translated_mesh(nested, s);
inner = carve(16, c2, 0.2, true);
outer = carve(16, c1, 0.25, true);
for k = 1:2
    for basis = ["dubiner", "peers"]
        nb = 2*nchoosek(k + 2, 2);
        if basis == "peers", nb = 4*nchoosek(k + 2, 2) + k + 1; end
        % nested: exact integral and contraction
        U = randn(nb, size(coarse.elements, 2));
        [V, info] = translated_field_between(sp, coarse, nested, k, k, U, s, basis = basis);
        I = sum(integral_of(nested, k, V, basis), 2);
        parent = containing_element(coarse, centroids(nested) - s);
        I_exact = sum(integral_of(nested, k, U(:, parent), basis, coarse, parent, s), 2);
        dI = max(abs(I - I_exact))/max(abs(I_exact));
        ratio = l2_norm(nested, k, V, basis)/l2_norm(coarse, k, U, basis);
        fprintf("d/e. nested   k = %d %-7s: integral %.1e, uncovered max %.1e, ||Pu||/||u|| = %.4f\n", ...
            k, basis, dI, max(info.uncovered), ratio);
        n_fail = n_fail + check(dI <= tol && max(info.uncovered) <= 1e-12 && ratio <= 1 + tol, ...
            sprintf("d/e. nested k = %d %s", k, basis));
        % snapped carves, inner inside outer + (c2 - c1)
        U = randn(nb, size(outer.elements, 2));
        [V, info] = translated_field_between(sp, outer, inner, k, k, U, c2 - c1, basis = basis);
        ratio = l2_norm(inner, k, V, basis)/l2_norm(outer, k, U, basis);
        fprintf("e.   snapped  k = %d %-7s: uncovered max %.1e, ||Pu||/||u|| = %.4f\n", ...
            k, basis, max(info.uncovered), ratio);
        n_fail = n_fail + check(max(info.uncovered) <= 1e-12 && ratio <= 1 + tol, ...
            sprintf("e. snapped k = %d %s", k, basis));
    end
end

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_shifted_projection_between:failed", "%d check(s) failed", n_fail);
end

%% local functions
function mesh = carve(N, centre, radius, snap)
% the solid side of the crossed background grid N x N split by a circle
curve = CircularCurve(centre, radius);
[~, mesh] = StructuredSimplicialGrid(N, N, "crossed").split(curve, max_normal_deviation = 90);
if snap
    mesh = snap_to_curve(mesh, curve, normal_sign = -1);
end
end

function mesh = translated_mesh(mesh, s)
mesh.coordinates(1:2, :) = mesh.coordinates(1:2, :) + s;
end

function d = rel(A, B)
d = max(abs(A - B), [], "all")/max(abs(B), [], "all");
end

function p = random_polynomial(k, m)
% m random polynomials of degree k: p(x) is n x m for x 2 x n
[i, j] = meshgrid(0:k);
e = [i(i + j <= k), j(i + j <= k)];
C = randn(size(e, 1), m);
ex = e(:, 1).';
ey = e(:, 2).';
p = @(x) ((x(1, :).').^ex .* (x(2, :).').^ey)*C;
end

function U = dubiner_dofs(mesh, k, p, s)
% the element-wise L2 projection onto Dubiner P_k of x -> p(x - s),
% component blocks
[rs, w] = get_2d_reference_quadrature(k + 2);
P = dubiner_basis(k, rs);
fit = (P.'*(w.*P)) \ (P.'.*w.');
[B, x3] = get_affine_map(mesh);
N = size(mesh.elements, 2);
U = [];
for K = 1:N
    v = p(B(:, :, K)*rs + x3(:, K) - s);
    U(:, K) = reshape(fit*v, [], 1); %#ok<AGROW>
end
end

function U = peers_dofs(mesh, k, q, s)
% PEERS dofs of the P_k tensor x -> q(x - s) (columns 11, 21, 12, 22)
U = dubiner_dofs(mesh, k, q, s);
U = [U; zeros(k + 1, size(U, 2))];
end

function values = field_values(k, U, r, G, basis)
% the field with dofs U (one column) at the reference points r of one
% element with inverse map G: n x m (Dubiner components) or n x 4 (PEERS)
switch basis
    case "dubiner"
        P = dubiner_basis(k, r);
        values = P*reshape(U, size(P, 2), []);
    case "peers"
        Psi = peers_points(k, r, repmat(G, 1, 1, size(r, 2)));
        values = reshape(pagemtimes(Psi, U), size(r, 2), 4);
end
end

function I = integral_of(mesh, k, V, basis, source_mesh, source, s)
% the integral over each element of mesh of the field V (columns per
% element); with a source mesh, the field of the source elements of
% source_mesh at x - s
[rs, w] = get_2d_reference_quadrature(k + 1);
[B, x3] = get_affine_map(mesh);
G = get_affine_map(mesh, inverse = true);
if nargin > 4
    [~, x3s] = get_affine_map(source_mesh);
    Gs = get_affine_map(source_mesh, inverse = true);
end
I = [];
for K = 1:size(mesh.elements, 2)
    area = abs(det(B(:, :, K)))/2;
    r = rs;  GK = G(:, :, K);
    if nargin > 4
        L = source(K);
        r = Gs(:, :, L)*(B(:, :, K)*rs + x3(:, K) - s - x3s(:, L));
        GK = Gs(:, :, L);
    end
    I(:, K) = area*(w.'*field_values(k, V(:, K), r, GK, basis)).'; %#ok<AGROW>
end
end

function n = l2_norm(mesh, k, V, basis)
[rs, w] = get_2d_reference_quadrature(k + 2);
B = get_affine_map(mesh);
G = get_affine_map(mesh, inverse = true);
n = 0;
for K = 1:size(mesh.elements, 2)
    v = field_values(k, V(:, K), rs, G(:, :, K), basis);
    n = n + abs(det(B(:, :, K)))/2*sum(w.*sum(v.^2, 2));
end
n = sqrt(n);
end

function c = centroids(mesh)
c = reshape(mean(reshape(mesh.coordinates(1:2, mesh.elements), 2, 3, []), 2), 2, []);
end

function el = containing_element(mesh, y)
% the element containing each point (interior points expected)
G = get_affine_map(mesh, inverse = true);
[~, x3] = get_affine_map(mesh);
el = zeros(1, size(y, 2));
for L = 1:size(mesh.elements, 2)
    l = G(:, :, L)*(y - x3(:, L));
    el(all(l > 1e-10, 1) & sum(l, 1) < 1 - 1e-10) = L;
end
assert(all(el > 0), "containing_element: a point in no element");
end

function [V, I] = composite_projection(old, new, kf, kt, U, s, basis, nsub, source)
% composite-quadrature projection: the new element split into nsub^2
% subtriangles, each with get_2d_reference_quadrature(max(kf, kt) + 2);
% x - s located in the old mesh, else the extrapolation source; also the
% integral of the shifted field over each new element
[rs0, w0] = get_2d_reference_quadrature(max(kf, kt) + 2);
[rs, w] = composite_rule(rs0, w0, nsub);
nP = numel(w);
[Bn, x3n] = get_affine_map(new);
Gn = get_affine_map(new, inverse = true);
[~, x3o] = get_affine_map(old);
Go = get_affine_map(old, inverse = true);
To = reshape(old.coordinates(1:2, old.elements), 2, 3, []);
lo = reshape(min(To, [], 2), 2, []);
hi = reshape(max(To, [], 2), 2, []);
if basis == "dubiner"
    Pt = dubiner_basis(kt, rs);
    fit = (Pt.'*(w.*Pt)) \ (Pt.'.*w.');
end
N = size(new.elements, 2);
V = [];  I = [];
for K = 1:N
    y = Bn(:, :, K)*rs + x3n(:, K) - s;
    el = zeros(1, nP);
    r = zeros(2, nP);
    candidates = find(all(lo <= max(y, [], 2) + 1e-12 & hi >= min(y, [], 2) - 1e-12, 1));
    for L = candidates
        open = find(el == 0);
        l = Go(:, :, L)*(y(:, open) - x3o(:, L));
        in = l(1, :) >= -1e-12 & l(2, :) >= -1e-12 & l(1, :) + l(2, :) <= 1 + 1e-12;
        el(open(in)) = L;
        r(:, open(in)) = l(:, in);
    end
    out = el == 0;
    el(out) = source(K);
    r(:, out) = Go(:, :, source(K))*(y(:, out) - x3o(:, source(K)));
    values = [];
    for L = unique(el)
        at = el == L;
        values(at, :) = field_values(kf, U(:, L), r(:, at), Go(:, :, L), basis); %#ok<AGROW>
    end
    I(:, K) = abs(det(Bn(:, :, K)))/2*(w.'*values).'; %#ok<AGROW>
    if basis == "dubiner"
        V(:, K) = reshape(fit*values, [], 1); %#ok<AGROW>
    else
        Phi = peers_points(kt, rs, repmat(Gn(:, :, K), 1, 1, nP));
        A = reshape(permute(Phi, [1 3 2]), nP*4, []);
        W = repmat(w, 4, 1);
        M = A.'*(W.*A);
        D = 1./sqrt(diag(M));
        V(:, K) = D.*((D.*M.*D.') \ (D.*(A.'*(W.*values(:))))); %#ok<AGROW>
    end
end
end

function [rs, w] = composite_rule(rs0, w0, n)
% the reference rule on each of the n^2 uniform subtriangles (weights sum 1)
lambda = [1 - rs0(1, :) - rs0(2, :); rs0(1, :); rs0(2, :)];
rs = zeros(2, 0);
w = zeros(0, 1);
for i = 0:n - 1
    for j = 0:n - 1 - i
        triangles = {[i j; i + 1 j; i j + 1]/n};
        if i + j <= n - 2
            triangles{end + 1} = [i + 1 j; i j + 1; i + 1 j + 1]/n; %#ok<AGROW>
        end
        for t = 1:numel(triangles)
            rs = [rs, triangles{t}.'*lambda]; %#ok<AGROW>
            w = [w; w0(:)/n^2]; %#ok<AGROW>
        end
    end
end
end

function state = random_state(mesh, k)
% a random Newmark solid state on mesh
d2 = nchoosek(k + 2, 2);
N = size(mesh.elements, 2);
F = size(mesh.faces, 2);
s = struct();
for level = ["", "_velocity", "_acceleration"]
    volume = ["displacement", "velocity", "acceleration"];
    volume = volume(["", "_velocity", "_acceleration"] == level);
    s.(volume) = randn(2*d2, N);
    s.(volume + "_skeleton") = randn(2*(k + 1), F);
    s.("stress" + level) = randn(4*d2 + k + 1, N);
    s.("rotation" + level) = randn(d2, N);
end
state = struct("solid", s);
end

function P = peers_points(k, r, G)
% PEERS_basis at the reference points r (2 x n), each with its own inverse
% affine map G (2 x 2 x n): n x nb x 4 (the 2 x 2 values column-major);
% the evaluation of LinearElasticityProblem's private shifted_projection
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
