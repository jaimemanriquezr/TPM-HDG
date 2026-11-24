% test_transfer_bc_darcy.m
% Transfer Path Method (Cockburn & Solano) boundary conditions for Darcy.
%
% Stage 1 (consistency): fitted unit square, identity mapping -> transfer
%   reduces to weak Dirichlet; solution must match the weak-BC solve.
% Stage 2 (unfitted): disk of radius R = 0.75 immersed in a background
%   grid of [-1,1]^2; cut elements discarded, Dirichlet data transferred
%   radially from the circle to the polygonal computational boundary.
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

dx = Measure(preset="dx");

%% ------------------------------------------------------------------ %%
%% Stage 1: fitted square, identity mapping == weak Dirichlet
%% ------------------------------------------------------------------ %%
mesh = StructuredSimplicialGrid(8, 8);

problem = FEProblem(Spaces= [FESpace([VPk, Pk], Condense=true), FESpace(Mk)]);
[lhs_op, flux_op] = DarcyMatrices(k, mesh, kappa);
problem.Operator = AffineOperator(lhs_op, flux_op);
[q, u, u_hat] = get_trial_functions(problem.Space); %#ok<ASGLU>
[r, w, v_hat] = get_test_functions(problem.Space);  %#ok<ASGLU>
F = assemble(f_ * w * dx, mesh);
problem.Operator = problem.Operator + F;

problem.BoundaryConditions = WeakEssentialBoundaryConditions(u_, "Dirichlet");
sol_weak = solve(problem, mesh, k);

problem.BoundaryConditions = TransferEssentialBoundaryConditions.darcy( ...
    u_, "Dirichlet", @(X) X, kappa);
sol_transfer = solve(problem, mesh, k);

diff_u = max(abs(sol_weak.Unknowns.field(2).dof - sol_transfer.Unknowns.field(2).dof), [], 'all');
fprintf("Stage 1 (identity mapping vs weak BC): max dof diff = %.3e\n", diff_u);
if diff_u > 1e-10
    fprintf("Stage 1 FAILED\n"); n_fail = n_fail + 1;
end

%% ------------------------------------------------------------------ %%
%% Stage 2: unfitted disk, radial transfer
%% ------------------------------------------------------------------ %%
R = 0.75;
phi = @(X) sqrt(X(1, :, :).^2 + X(2, :, :).^2) - R;
radial_map = @(X) R * X ./ vecnorm(X);

n_ref = 4;
[err_u, err_q, h_list] = deal(zeros(n_ref, 1));
for i = 1:n_ref
    n = 8 * 2^(i-1);
    h_list(i) = 2/n;
    background = StructuredSimplicialGrid(n, n, -1, 1, -1, 1);
    mesh = immerse_mesh(background, phi);

    problem = FEProblem(Spaces= [FESpace([VPk, Pk], Condense=true), FESpace(Mk)]);
    [lhs_op, flux_op] = DarcyMatrices(k, mesh, kappa);
problem.Operator = AffineOperator(lhs_op, flux_op);
    F = assemble(f_ * w * dx, mesh);
    problem.Operator = problem.Operator + F;
    problem.BoundaryConditions = TransferEssentialBoundaryConditions.darcy( ...
        u_, "Dirichlet", radial_map, kappa);

    solution = solve(problem, mesh, k);
    u_h = solution.Unknowns.field(2);
    q_h = solution.Unknowns.field(1);
    err_u(i) = compute_l2_error(u_h, mesh, u_);
    err_q(i) = compute_l2_error(q_h, mesh, Datum(q_sym));
    fprintf("n = %3d, h = %.4f, L2 err u = %.4e, L2 err q = %.4e\n", ...
        n, h_list(i), err_u(i), err_q(i));
end
eoc_u = compute_eoc(err_u, h_list);
eoc_q = compute_eoc(err_q, h_list);
fprintf("EOC(u): %s (expected ~ %d)\n", num2str(eoc_u.', "%.2f "), k+1);
fprintf("EOC(q): %s (expected ~ %d)\n", num2str(eoc_q.', "%.2f "), k+1);

% rates between the coarsest and finest levels (the mean EOC): the
% per-level EOCs of the carved disk oscillate around k+1
rate_u = log(err_u(1)/err_u(end)) / log(h_list(1)/h_list(end));
rate_q = log(err_q(1)/err_q(end)) / log(h_list(1)/h_list(end));
fprintf("rate(u), rate(q) over all levels: %.2f, %.2f (pass >= %.2f)\n", ...
    rate_u, rate_q, k + 1 - 0.25);
if min(rate_u, rate_q) < k + 1 - 0.25
    fprintf("Stage 2 rate FAILED\n"); n_fail = n_fail + 1;
end

fprintf("\n%d check(s) failed.\n", n_fail);
