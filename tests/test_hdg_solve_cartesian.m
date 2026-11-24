%% hdg_solve_cartesian: Q^k HDG for advection-diffusion on a CartesianGrid
% With beta != 0 and tau = epsilon + |beta.n|, a u in Q^2 is reproduced
% exactly on a non-uniform grid, and for the data of the LowRankHDG note
% (epsilon = 1, beta = (1, 0.5), u = sin(pi x) sin(pi y) + exp(x y)/2)
% the L2 error of u_h decreases at rate k + 1 (k = 1:3). With a reaction
% term the patch test still holds, mass_rhs agrees with the same source
% given as a function, and implicit Euler steps (reaction = 1/dt,
% mass_rhs = u_old/dt) converge at first order in dt.
clc; clear;
n_fail = 0;
syms x y real

mesh = CartesianGrid([0 0.1 0.35 0.5 1], [-1 -0.2 0.3 2]);
[u, f] = manufactured(x^2*y - 2*x*y^2 + x - y, 0.3, [1; -0.5], x, y);
s = hdg_solve_cartesian(mesh, 2, f, u, epsilon = 0.3, beta = [1; -0.5]);
E = QkBasis(2).compute_l2_error(s.u, mesh, exact_solution = u);
n_fail = n_fail + check("advective patch test", E, 1e-12);

[u, f] = manufactured(sin(pi*x)*sin(pi*y) + exp(x*y)/2, 1, [1; 0.5], x, y);
for k = 1:3
    E = zeros(1, 3);
    for n = 1:3
        g = CartesianGrid(4*2^n);
        s = hdg_solve_cartesian(g, k, f, u, epsilon = 1, beta = [1; 0.5]);
        E(n) = QkBasis(k).compute_l2_error(s.u, g, exact_solution = u);
    end
    fprintf("  k = %d errors %s\n", k, mat2str(E, 3));
    n_fail = n_fail + check(sprintf("rate k = %d", k), abs(log2(E(2)/E(3)) - (k + 1)), 0.15);
end
% reaction: Q^2 patch test with c = 2.5
[u, f] = manufactured(x^2*y - 2*x*y^2 + x - y, 0.3, [1; -0.5], x, y, 2.5);
s = hdg_solve_cartesian(mesh, 2, f, u, epsilon = 0.3, beta = [1; -0.5], reaction = 2.5);
E = QkBasis(2).compute_l2_error(s.u, mesh, exact_solution = u);
n_fail = n_fail + check("reaction patch test", E, 1e-12);

% mass_rhs: f + w with w in Q^2 given by coefficients equals the solve
% with the source f + w as a function handle
w_fun = u; w_coef = s.u;                       % s.u reproduces u exactly
s1 = hdg_solve_cartesian(mesh, 2, @(x, y) f(x, y) + w_fun(x, y), u, epsilon = 0.3, ...
    beta = [1; -0.5], reaction = 2.5);
s2 = hdg_solve_cartesian(mesh, 2, f, u, epsilon = 0.3, beta = [1; -0.5], ...
    reaction = 2.5, mass_rhs = w_coef);
n_fail = n_fail + check("mass_rhs = source", norm(s1.u - s2.u, 'fro') / norm(s1.u, 'fro'), 1e-12);

% implicit Euler for u_t - lap(u) + beta.grad(u) = f to T = 0.2 (k = 3,
% N = 16, spatial error negligible): first order in dt
syms t real
u_st = exp(-t)*(sin(pi*x)*sin(pi*y) + exp(x*y)/2);
f_st = diff(u_st, t) - (diff(u_st, x, 2) + diff(u_st, y, 2)) + diff(u_st, x) + 0.5*diff(u_st, y);
u_t = matlabFunction(u_st, Vars = [t x y]); f_t = matlabFunction(f_st, Vars = [t x y]);
g16 = CartesianGrid(16); k = 3; T_end = 0.2;
[xr, wr] = get_2d_reference_quadrature(k + 1, "Quad");
P0 = QkBasis.evaluate_basis(k, xr);
xq = g16.get_physical_points(xr);
E = zeros(1, 3);
for i = 1:3
    dt = 0.05 / 2^(i - 1);
    u_old = P0.' * (wr .* reshape(u_t(0, xq(1, :, :), xq(2, :, :)), [], size(xq, 3)));
    for n = 1:round(T_end/dt)
        tn = n*dt;
        s = hdg_solve_cartesian(g16, k, @(x, y) f_t(tn, x, y) + 0*x, @(x, y) u_t(tn, x, y) + 0*x, ...
            beta = [1; 0.5], reaction = 1/dt, mass_rhs = u_old/dt);
        u_old = s.u;
    end
    E(i) = QkBasis(k).compute_l2_error(s.u, g16, exact_solution = @(x, y) u_t(T_end, x, y) + 0*x);
end
fprintf("  implicit Euler errors %s\n", mat2str(E, 3));
n_fail = n_fail + check("implicit Euler rate", abs(log2(E(2)/E(3)) - 1), 0.1);
fprintf("%d check(s) failed\n", n_fail);

function [u, f] = manufactured(u_sym, epsilon, beta, x, y, c)
if nargin < 6
    c = 0;
end
f_sym = -epsilon*(diff(u_sym, x, 2) + diff(u_sym, y, 2)) ...
        + beta(1)*diff(u_sym, x) + beta(2)*diff(u_sym, y) + c*u_sym;
u0 = matlabFunction(u_sym, Vars = [x y]);
f0 = matlabFunction(f_sym, Vars = [x y]);
u = @(x, y) u0(x, y) + 0*x;
f = @(x, y) f0(x, y) + 0*x;
end

function failed = check(name, err, tol)
    failed = ~(err <= tol);
    fprintf("  %-22s err = %.2e  %s\n", name, err, string(~failed));
end
