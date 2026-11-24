% test_weak_bc_baseline.m
% Baseline for the declarative solve path: fitted Darcy on the unit square
% with weakly imposed Dirichlet data, DarcyMatrices local problem.
clc; clear;
n_fail = 0;

kappa = eye(2);
k = 2;

syms x y
u_sym = sin(pi*x) .* sin(pi*y);
q_sym = -kappa * gradient(u_sym, [x y]);
f_sym = divergence(q_sym, [x y]);

u_ = Datum(u_sym);
f_ = Datum(f_sym);

Pk = DubinerBasis(k);
VPk = Pk.to_vector();
Mk = LegendreBasis(k);

problem = FEProblem(Spaces= [FESpace([VPk, Pk], Condense=true), FESpace(Mk)]);
problem.BoundaryConditions = WeakEssentialBoundaryConditions(u_, "Dirichlet");

[q, u, u_hat] = get_trial_functions(problem.Space);
[r, w, v_hat] = get_test_functions(problem.Space);
dx = Measure(preset="dx");

n_ref = 3;
[err_u, h_list] = deal(zeros(n_ref, 1));
for i = 1:n_ref
    n = 4 * 2^(i-1);
    h_list(i) = 1/n;
    mesh = StructuredSimplicialGrid(n, n);

    problem.Mesh = mesh;
    F = assemble(f_ * w * dx, mesh);
    [lhs_op, flux_op] = DarcyMatrices(k, mesh, kappa);
problem.Operator = AffineOperator(lhs_op, flux_op);
    problem.Operator = problem.Operator + F;

    solution = solve(problem, mesh, k);
    u_h = solution.Unknowns.field(2);
    err_u(i) = compute_l2_error(u_h, mesh, u_);
    fprintf("n = %2d, L2 err u = %.4e\n", n, err_u(i));
end
eoc = compute_eoc(err_u, h_list);
fprintf("EOC(u): %s (expected ~ %d)\n", num2str(eoc.', "%.2f "), k+1);

% rate between the coarsest and finest levels (the mean EOC)
rate = log(err_u(1)/err_u(end)) / log(h_list(1)/h_list(end));
fprintf("rate(u) over all levels: %.2f (pass >= %.2f)\n", rate, k + 1 - 0.25);
if rate < k + 1 - 0.25
    fprintf("rate(u) FAILED\n"); n_fail = n_fail + 1;
end

fprintf("\n%d check(s) failed.\n", n_fail);
