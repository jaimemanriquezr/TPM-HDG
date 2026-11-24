% test_translated_field_extension.m
% The option extension_degree = q of
% LinearElasticityProblem.translated_field_between: on the part R_K of a
% new element K that the old mesh does not cover, the projection
% integrates the degree-q Taylor polynomial of the source element's
% polynomial about x_K, the point of the source element nearest to the
% centroid of R_K. Meshes: the snapped carves of the disc of radius 0.2
% (crossed grid, N = 16) with centres (0.5, 0.5) and (0.53, 0.5).
%   (a) gate: no option and extension_degree = [] give the same outputs
%       (bit for bit);
%   (b) exactness: with q = 1 a global polynomial of degree <= 1, and with
%       q = 0 a constant, are reproduced (1e-12 relative) from either carve
%       onto the other, k = 1, 2, 3; max(info.uncovered) > 0.05;
%   (c) for a global polynomial p of degree k = 2, 3 and q = 1, on every
%       element with uncovered fraction > 1e-6: the result differs from the
%       full extrapolation (which reproduces p) and equals (1e-11) the L2
%       projection of {p on the covered part, Taylor_1(p; x_K) on R_K},
%       computed independently (R_K by polyshape boolean operations, its
%       triangulation and a Grundmann-Moeller rule; x_K from the
%       polyshape centroid of R_K, 1e-12);
%   (d) q = 1, a field of degree 1, shift (0.01, 0.005): exact (1e-12).
failures = 0;
N = 16; radius = 0.2;
carve = @(c) carve_solid(local_background(N, c, radius), CircularCurve(c, radius));
meshes = {carve([0.5; 0.5]), carve([0.53; 0.5])};

% (a) gate
k = 2;
rng(2); U = randn(4*nchoosek(k + 2, 2), size(meshes{1}.elements, 2));
sp = LinearElasticityProblem(k);
[V0, info0] = sp.translated_field_between(meshes{1}, meshes{2}, k, k, U, [0.01; 0.005]);
[V1, info1] = sp.translated_field_between(meshes{1}, meshes{2}, k, k, U, [0.01; 0.005], extension_degree = []);
fprintf("(a) no option vs extension_degree = []: isequal V %d, info %d\n", isequal(V0, V1), isequal(info0, info1));
failures = failures + check(isequal(V0, V1) && isequal(info0, info1), "(a) gate");

% (b) exactness for affine fields (q = 1) and constants (q = 0)
affine = {@(x, y) 0.3 + 1.2*x - 0.7*y, @(x, y) -0.4 + 0.5*x + 0.9*y};
constant = {@(x, y) 0.7 + 0*x, @(x, y) -1.1 + 0*x};
for k = 1:3
    sp = LinearElasticityProblem(k);
    for pair = [1 2; 2 1]
        old = meshes{pair(1)}; new = meshes{pair(2)};
        for q = [1 0]
            if q == 1, f = affine; else, f = constant; end
            [V, info] = sp.translated_field_between(old, new, k, k, L2proj2D(old, f, k), [0; 0], ...
                extension_degree = q);
            ref = L2proj2D(new, f, k);
            err = max(abs(V - ref), [], "all")/max(abs(ref), [], "all");
            fprintf("(b) k = %d, carve %d -> %d, q = %d: max uncovered %.3f, relative error %.2e\n", ...
                k, pair, q, max(info.uncovered), err);
            failures = failures + check(err <= 1e-12 && max(info.uncovered) > 0.05, ...
                sprintf("(b) k = %d, %d -> %d, q = %d", k, pair, q));
        end
    end
end

% (c) the extension is the claimed one
old = meshes{1}; new = meshes{2};
T_old = reshape(old.coordinates(1:2, old.elements), 2, 3, []);
covered = union(arrayfun(@(i) polyshape(T_old(1, :, i), T_old(2, :, i)), 1:size(T_old, 3)));
[B, x3] = get_affine_map(new);
for k = 2:3
    sp = LinearElasticityProblem(k);
    [p, grad_p] = random_polynomial(k, 20 + k);
    f = {p, @(x, y) p(y, x)};
    U = L2proj2D(old, f, k);
    V_full = sp.translated_field_between(old, new, k, k, U, [0; 0]);
    [V, info] = sp.translated_field_between(old, new, k, k, U, [0; 0], extension_degree = 1);
    [rs, w] = get_2d_reference_quadrature(k + 2);
    d2 = nchoosek(k + 2, 2);
    elements = find(info.uncovered > 1e-6);
    [err, err_point, change] = deal(0);
    for K = elements
        X = new.coordinates(1:2, new.elements(:, K));
        R = subtract(polyshape(X(1, :), X(2, :)), covered);
        [cx, cy] = centroid(R);
        E = old.coordinates(1:2, old.elements(:, info.source(K)));
        x_K = nearest_in_triangle(E, [cx; cy]);
        err_point = max(err_point, norm(x_K - info.taylor_point(:, K)));
        % moments: (p, w)_K - (p - Taylor_1 p, w)_R
        basis = @(x) dubiner_basis(k, B(:, :, K) \ (x - x3(:, K)));
        [xq, wq] = triangle_rule(X, rs, w);
        Phi = basis(xq);
        mass = Phi.'*(wq.*Phi);
        m = Phi.'*(wq.*[f{1}(xq(1, :), xq(2, :)).', f{2}(xq(1, :), xq(2, :)).']);
        tri = triangulation(R);
        for t = 1:size(tri.ConnectivityList, 1)
            [xq, wq] = triangle_rule(tri.Points(tri.ConnectivityList(t, :), :).', rs, w);
            g1 = grad_p(x_K(1), x_K(2));
            g2 = flip(grad_p(x_K(2), x_K(1)));
            rest = [p(xq(1, :), xq(2, :)).' - (p(x_K(1), x_K(2)) + (xq - x_K).'*g1), ...
                    p(xq(2, :), xq(1, :)).' - (p(x_K(2), x_K(1)) + (xq - x_K).'*g2)];
            m = m - basis(xq).'*(wq.*rest);
        end
        ref = mass \ m;
        V_K = reshape(V(:, K), d2, 2);
        err = max(err, max(abs(V_K - ref), [], "all")/max(abs(ref), [], "all"));
        change = max(change, max(abs(V(:, K) - V_full(:, K)))/max(abs(V_full(:, K))));
    end
    fprintf("(c) k = %d, %d elements: |x_K - polyshape x_K| = %.2e, |V - brute force|/|V| = %.2e, " + ...
        "|V - V_full|/|V_full| = %.2e\n", k, numel(elements), err_point, err, change);
    failures = failures + check(err_point <= 1e-12, sprintf("(c) k = %d: Taylor point", k));
    failures = failures + check(err <= 1e-11, sprintf("(c) k = %d: brute-force projection", k));
    failures = failures + check(change > 1e-6, sprintf("(c) k = %d: differs from the full extrapolation", k));
end

% (d) shift nonzero
shift = [0.01; 0.005];
for k = 1:3
    sp = LinearElasticityProblem(k);
    [V, info] = sp.translated_field_between(meshes{1}, meshes{2}, k, k, L2proj2D(meshes{1}, affine, k), ...
        shift, extension_degree = 1);
    ref = L2proj2D(meshes{2}, cellfun(@(g) @(x, y) g(x - shift(1), y - shift(2)), affine, ...
        "UniformOutput", false), k);
    err = max(abs(V - ref), [], "all")/max(abs(ref), [], "all");
    fprintf("(d) k = %d, shift (0.01, 0.005): max uncovered %.3f, relative error %.2e\n", ...
        k, max(info.uncovered), err);
    failures = failures + check(err <= 1e-12, sprintf("(d) k = %d", k));
end

if failures == 0
    fprintf("test_translated_field_extension: all checks passed\n");
else
    fprintf("test_translated_field_extension: %d check(s) failed\n", failures);
end

function grid = local_background(N, center, radius)
% the cells of StructuredSimplicialGrid(N, N, "crossed") near the disc
lo = max(floor((center - radius)*N) - 2, 0);
hi = min(ceil((center + radius)*N) + 2, N);
grid = StructuredSimplicialGrid(hi(1) - lo(1), hi(2) - lo(2), lo(1)/N, hi(1)/N, lo(2)/N, hi(2)/N, ...
    style = "crossed");
end

function mesh = carve_solid(background, curve)
% the inside of the circle, cut elements discarded, interface vertices snapped
[~, mesh] = background.split(curve, max_normal_deviation = 90);
mesh = snap_to_curve(mesh, curve, normal_sign = -1);
end

function [p, grad_p] = random_polynomial(k, seed)
% a polynomial of total degree k with random coefficients, and its gradient
[a, b] = meshgrid(0:k);
keep = a + b <= k;
a = a(keep).'; b = b(keep).';
rng(seed); c = randn(1, numel(a));
p = @(x, y) reshape(sum(c.*x(:).^a.*y(:).^b, 2), size(x));
grad_p = @(x, y) [sum(c.*a.*x.^max(a - 1, 0).*y.^b); sum(c.*b.*x.^a.*y.^max(b - 1, 0))];
end

function [x, w] = triangle_rule(X, rs, w)
% the reference rule (sum(w) = 1) on the triangle with vertices X (2 x 3)
x = X(:, 1) + (X(:, 2) - X(:, 1)).*rs(1, :) + (X(:, 3) - X(:, 1)).*rs(2, :);
area = abs(det([X(:, 2) - X(:, 1), X(:, 3) - X(:, 1)]))/2;
w = area*w(:);
end

function x = nearest_in_triangle(X, c)
% the point of the closed triangle X (2 x 3) nearest to c
if isinterior(polyshape(X(1, :), X(2, :)), c(1), c(2))
    x = c;
    return
end
best = inf;
for e = 1:3
    a = X(:, e); d = X(:, mod(e, 3) + 1) - a;
    y = a + min(max((c - a).'*d/(d.'*d), 0), 1)*d;
    if norm(c - y) < best, best = norm(c - y); x = y; end
end
end
