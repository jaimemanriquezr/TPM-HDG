%% QkBasis on CartesianGrid
% On a non-uniform grid, for k = 0:4: the element mass matrix is
% volume * I (orthonormal basis), the face mass matrix of LegendreBasis(k)
% is area * I; the L2 projection of u in Q^k is exact (compute_l2_error,
% gradients via evaluate_gradient, values at physical points); face and
% element bases pass the compatibility checks, a Dubiner basis does not.
% On uniform grids, the L2-projection error of a smooth function
% decreases at rate k + 1.
clc; clear;
n_fail = 0;
tol = 1e-11;

mesh = CartesianGrid([0 0.1 0.35 0.5 1], [-1 -0.2 0.3 2]);
n_el = size(mesh.elements, 2);
for k = 0:4
    basis = QkBasis(k);
    d = basis.Dimension;
    fprintf("k = %d\n", k);

    M = mass_matrix(mesh, basis);
    n_fail = n_fail + check("element mass", ...
        max(abs(M - reshape(mesh.volume, 1, 1, []) .* eye(d)), [], 'all'), tol);
    Mf = mass_matrix(mesh, LegendreBasis(k));
    n_fail = n_fail + check("face mass", ...
        max(abs(Mf - reshape(mesh.area, 1, 1, []) .* eye(k + 1)), [], 'all'), tol);

    % u in Q^k and its gradient
    c = double(k > 0);                               % the linear part is in Q^k for k >= 1
    u = @(x, y) (x - 0.3).^k .* (y + 0.5).^k + c*(x - 2*y);
    ux = @(x, y) k*(x - 0.3).^max(k-1, 0) .* (y + 0.5).^k + c;
    uy = @(x, y) k*(x - 0.3).^k .* (y + 0.5).^max(k-1, 0) - 2*c;
    dof = project(basis, mesh, u);
    E = basis.compute_l2_error(dof, mesh, exact_solution = u);
    n_fail = n_fail + check("Q^k projection", E, tol);

    rs = [0.13 0.5 0.92; 0.71 0.05 0.5];
    X = mesh.get_physical_points(rs);
    dp = basis.evaluate_gradient(rs, mesh);          % n_q x d x 2 x 1 x n_el
    gx = squeeze(pagemtimes(dp(:, :, 1, 1, :), reshape(dof, d, 1, 1, 1, [])));
    gy = squeeze(pagemtimes(dp(:, :, 2, 1, :), reshape(dof, d, 1, 1, 1, [])));
    ex = squeeze(ux(X(1, :, :), X(2, :, :)));
    ey = squeeze(uy(X(1, :, :), X(2, :, :)));
    n_fail = n_fail + check("gradient", max(abs([gx - ex; gy - ey]), [], 'all'), 1e-9);

    % values at physical points (element search)
    P = reshape(X, 2, []);
    [~, el] = mesh.get_2d_local(P);
    vals = basis.evaluate(P, mesh, coordinate_system = "physical");
    uh = sum(vals .* dof(:, el).', 2);
    n_fail = n_fail + check("physical eval", max(abs(uh - u(P(1, :), P(2, :)).')), tol);
end

% compatibility rules for face bases
ok_face = true; ok_simplex = false;
try, mustBeSameElement(LegendreBasis(2), mesh); mustBeCompatible(LegendreBasis(1), LegendreBasis(1)); catch, ok_face = false; end
try, mustBeSameElement(DubinerBasis(2), mesh); catch, ok_simplex = true; end
n_fail = n_fail + check("face basis allowed", double(~ok_face), 0);
n_fail = n_fail + check("Dubiner rejected", double(~ok_simplex), 0);
[s1, w1] = get_quadrature([QkBasis(2), LegendreBasis(2)], 1);
% order 2 -> get_quadrature asks for k + 1 = 3 -> 4 Gauss points
n_fail = n_fail + check("mixed face rule", abs(numel(w1) - 4) + abs(sum(w1) - 1) + (size(s1, 1) ~= 1), tol);

% convergence of the L2 projection
f = @(x, y) sin(pi*x) .* exp(y);
Ns = [4 8 16 32];
fprintf("rates of the L2 projection of sin(pi x) exp(y):\n");
for k = 1:4
    E = zeros(size(Ns));
    for n = 1:numel(Ns)
        g = CartesianGrid(Ns(n));
        E(n) = QkBasis(k).compute_l2_error(project(QkBasis(k), g, f), g, exact_solution = f);
    end
    rates = log2(E(1:end-1) ./ E(2:end));
    fprintf("  k = %d  errors %s  rates %s\n", k, mat2str(E, 3), mat2str(rates, 3));
    n_fail = n_fail + check(sprintf("rate k = %d", k), abs(rates(end) - (k + 1)), 0.2);
end
fprintf("%d check(s) failed\n", n_fail);

function dof = project(basis, mesh, u)
% L2 projection, d x n_el
[xy, rs, w] = quadpoints_elements_2d(mesh, basis.Order + 1, output="array");
P = basis.evaluate(rs);                               % n_q x d
U = reshape(u(xy(1, :, :), xy(2, :, :)), [], size(mesh.elements, 2));
b = mesh.volume .* (P.' * (w .* U));                  % d x n_el
M = mass_matrix(mesh, basis);
dof = reshape(pagemldivide(M, reshape(b, [], 1, size(b, 2))), [], size(b, 2));
end

function failed = check(name, err, tol)
    failed = ~(err <= tol);
    fprintf("  %-20s err = %.2e  %s\n", name, err, string(~failed));
end
