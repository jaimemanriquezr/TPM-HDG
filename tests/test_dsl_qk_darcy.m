%% Q^k Darcy on a CartesianGrid in the form language
% The HDG Darcy system of dsl_darcy_example, with QkBasis / LegendreBasis
% on a CartesianGrid, is solved on the declarative path and compared with
% the standalone hdg_poisson_cartesian (same formulation, tau = 1): the
% coefficients of u_h and q_h agree to machine precision for k = 1:3 on a
% uniform and on a non-uniform grid, and the L2 error of u_h on uniform
% grids decreases at rate k + 1.
clc; clear;
n_fail = 0;
tol = 1e-12;

ufun = @(x, y) sin(pi*x).*sin(pi*y);
ffun = @(x, y) 2*pi^2*sin(pi*x).*sin(pi*y);
g = Datum(ufun); f = Datum(ffun);

grids = {CartesianGrid(6), CartesianGrid([0 0.1 0.35 0.5 0.8 1], [0 0.3 0.45 1])};
names = ["uniform", "non-uniform"];
for gi = 1:numel(grids)
    mesh = grids{gi}.add_stabilization(1);
    for k = 1:3
        sol = solve_dsl(mesh, k, f, g);
        ref = hdg_poisson_cartesian(mesh, k, ffun, ufun);
        du = max(abs(sol.Unknowns.field(2).dof - ref.u), [], 'all');
        dq = max(abs(sol.Unknowns.field(1).dof - ref.q), [], 'all');
        n_fail = n_fail + check(sprintf("%s k = %d u", names(gi), k), du, tol);
        n_fail = n_fail + check(sprintf("%s k = %d q", names(gi), k), dq, tol);
    end
end

% convergence of u_h
Ns = [4 8 16];
for k = 1:3
    E = zeros(size(Ns));
    for n = 1:numel(Ns)
        mesh = CartesianGrid(Ns(n)).add_stabilization(1);
        sol = solve_dsl(mesh, k, f, g);
        E(n) = compute_l2_error(sol.Unknowns.field(2), mesh, g);
    end
    rates = log2(E(1:end-1) ./ E(2:end));
    fprintf("  k = %d  errors %s  rates %s\n", k, mat2str(E, 3), mat2str(rates, 3));
    n_fail = n_fail + check(sprintf("rate k = %d", k), abs(rates(end) - (k + 1)), 0.15);
end
fprintf("%d check(s) failed\n", n_fail);

function sol = solve_dsl(mesh, k, f, g)
    Qk = QkBasis(k); VQk = Qk.to_vector(); Mk = LegendreBasis(k);
    problem = FEProblem(Spaces = [FESpace([VQk, Qk], Condense = true), FESpace(Mk)]);
    problem.BoundaryConditions = WeakEssentialBoundaryConditions(g, "Dirichlet");
    [q, u, u_hat] = get_trial_functions(problem.Space);
    [r, w, v_hat] = get_test_functions(problem.Space);
    dx = Measure(preset = "dx");
    dS = Measure(preset = "dS");
    n = ExteriorNormal();
    problem.Operator = AffineOperator( ...
          dot(q, r)*dx - u*divergence(r)*dx + u_hat*dot(r, n)*dS ...
        + divergence(q)*w*dx + u*w*dS - u_hat*w*dS ...
        + f*w*dx ...
        - (dot(q, n)*v_hat*dS + u*v_hat*dS - u_hat*v_hat*dS));
    sol = solve(problem, mesh, k);
end

function failed = check(name, err, tol)
    failed = ~(err <= tol);
    fprintf("  %-20s err = %.2e  %s\n", name, err, string(~failed));
end
