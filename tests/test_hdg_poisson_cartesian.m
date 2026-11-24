%% hdg_poisson_cartesian: Q^k HDG for Poisson on a CartesianGrid
% A harmonic u in Q^2 (not a product of 1-D polynomials) is reproduced
% exactly on a non-uniform grid; for u = sin(pi x) sin(pi y) the L2 error
% of u_h decreases at rate k + 1 (k = 0:3).
clc; clear;
n_fail = 0;

mesh = CartesianGrid([0 0.2 0.7 1], [0 0.4 1]);
up = @(x, y) x.^2 - y.^2 + x.*y;
s = hdg_poisson_cartesian(mesh, 2, @(x, y) 0*x, up);
E = QkBasis(2).compute_l2_error(s.u, mesh, exact_solution = up);
n_fail = n_fail + check("Q^2 patch test", E, 1e-12);

u = @(x, y) sin(pi*x).*sin(pi*y);
f = @(x, y) 2*pi^2*u(x, y);
for k = 0:3
    E = zeros(1, 3);
    for n = 1:3
        g = CartesianGrid(4*2^n);
        s = hdg_poisson_cartesian(g, k, f, u);
        E(n) = QkBasis(k).compute_l2_error(s.u, g, exact_solution = u);
    end
    rate = log2(E(end-1)/E(end));
    fprintf("  k = %d errors %s\n", k, mat2str(E, 3));
    n_fail = n_fail + check(sprintf("rate k = %d", k), abs(rate - (k + 1)), 0.15);
end
fprintf("%d check(s) failed\n", n_fail);

function failed = check(name, err, tol)
    failed = ~(err <= tol);
    fprintf("  %-16s err = %.2e  %s\n", name, err, string(~failed));
end
