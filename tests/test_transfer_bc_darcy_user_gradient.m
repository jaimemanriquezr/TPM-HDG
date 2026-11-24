% test_transfer_bc_darcy_user_gradient.m
% Transfer Path Method BCs for Darcy with a user-defined gradient
% reconstruction: no preset factory (.darcy) — the script declares the
% constitutive relation grad(u) = -kappa^{-1} q itself, in both forms the
% path_gradient option accepts: the symbolic expression -kappa \ q in the
% flux trial function, and a low-level function handle (local function at
% the bottom of the script).
%
% Stage 1 (consistency): fitted unit square, identity mapping -> transfer
%   reduces to weak Dirichlet, and must also match the .darcy factory
%   dof-for-dof.
% Stage 2 (unfitted): disk of radius R = 0.75 immersed in [-1,1]^2,
%   Dirichlet data transferred radially; EOC ~ k+1.
clc; clear;
n_fail = 0;

% Non-symmetric on purpose: kappa^{-T} ~= kappa^{-1}, so a transpose bug
% in the reconstruction degrades the stage-2 EOC. Unphysical as a
% permeability, but the PDE stays elliptic (symmetric part is PD).
kappa = [2 1; 0 3];
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

% The user-declared reconstruction: fix the PDE parameters by closure and
% match the documented path_gradient signature (see the class help).
pg = @(problem, mesh, rs, K_list, m, n_local) ...
    darcy_user_gradient(problem, rs, m, n_local, kappa);

%% ------------------------------------------------------------------ %%
%% Stage 1: fitted square — identity mapping == weak Dirichlet,
%%          and user gradient == .darcy factory
%% ------------------------------------------------------------------ %%
mesh = StructuredSimplicialGrid(8, 8);

problem = FEProblem(Spaces= [FESpace([VPk, Pk], Condense=true), FESpace(Mk)]);
[lhs_op, flux_op] = DarcyMatrices(k, mesh, kappa);
problem.Operator = AffineOperator(lhs_op, flux_op);
[q, u, ~] = get_trial_functions(problem.Space);
[~, w, ~] = get_test_functions(problem.Space);

% Fully declarative alternative: state the constitutive law itself as an
% expression in the flux trial function; the contraction with the
% transfer direction happens inside the BC class.
pg_expr = -kappa \ q;

% Multi-term sum (LinearCombination lowering): the two q-terms have
% different wrapper depths and add up to -kappa^{-1} q; the zero u-term
% exercises the second-field dof offset.
pg_sum = -((2*kappa) \ q) + q*(-0.5*inv(kappa)) + [0; 0]*u;
F = assemble(f_ * w * dx, mesh);
problem.Operator = problem.Operator + F;

problem.BoundaryConditions = WeakEssentialBoundaryConditions(u_, "Dirichlet");
sol_weak = solve(problem, mesh, k);

problem.BoundaryConditions = TransferEssentialBoundaryConditions( ...
    u_, "Dirichlet", @(X) X, path_gradient = pg);
sol_user = solve(problem, mesh, k);

problem.BoundaryConditions = TransferEssentialBoundaryConditions.darcy( ...
    u_, "Dirichlet", @(X) X, kappa);
sol_factory = solve(problem, mesh, k);

problem.BoundaryConditions = TransferEssentialBoundaryConditions( ...
    u_, "Dirichlet", @(X) X, path_gradient = pg_expr);
sol_expr = solve(problem, mesh, k);

problem.BoundaryConditions = TransferEssentialBoundaryConditions( ...
    u_, "Dirichlet", @(X) X, path_gradient = pg_sum);
sol_sum = solve(problem, mesh, k);

diff_weak    = max(abs(sol_weak.Unknowns.field(2).dof    - sol_user.Unknowns.field(2).dof), [], 'all');
diff_factory = max(abs(sol_factory.Unknowns.field(2).dof - sol_user.Unknowns.field(2).dof), [], 'all');
diff_expr    = max(abs(sol_expr.Unknowns.field(2).dof    - sol_user.Unknowns.field(2).dof), [], 'all');
diff_sum     = max(abs(sol_sum.Unknowns.field(2).dof     - sol_expr.Unknowns.field(2).dof), [], 'all');
fprintf("Stage 1 (identity mapping vs weak BC):   max dof diff = %.3e\n", diff_weak);
fprintf("Stage 1 (user gradient vs .darcy):       max dof diff = %.3e\n", diff_factory);
fprintf("Stage 1 (-kappa\\q expression vs handle): max dof diff = %.3e\n", diff_expr);
fprintf("Stage 1 (multi-term sum vs -kappa\\q):    max dof diff = %.3e\n", diff_sum);
if diff_weak > 1e-10 || diff_factory > 1e-12 || diff_expr > 1e-12 || diff_sum > 1e-12
    fprintf("Stage 1 FAILED\n"); n_fail = n_fail + 1;
end

%% ------------------------------------------------------------------ %%
%% Stage 2: unfitted disk, radial transfer with the user gradient
%% ------------------------------------------------------------------ %%
R = 0.75;
phi = @(X) sqrt(X(1, :, :).^2 + X(2, :, :).^2) - R;
radial_map = @(X) R * X ./ vecnorm(X);

n_ref = 3;
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
    problem.BoundaryConditions = TransferEssentialBoundaryConditions( ...
        u_, "Dirichlet", radial_map, path_gradient = pg_expr);

    solution = solve(problem, mesh, k);
    if i == 1
        % nonzero transfer paths: the real probe of the sum lowering
        % (stage 1 has m = 0, so the gradient never enters there)
        problem.BoundaryConditions = TransferEssentialBoundaryConditions( ...
            u_, "Dirichlet", radial_map, path_gradient = pg_sum);
        sol_sum2 = solve(problem, mesh, k);
        diff_sum2 = max(abs(sol_sum2.Unknowns.field(2).dof ...
                          - solution.Unknowns.field(2).dof), [], 'all');
        fprintf("Stage 2 (multi-term sum vs -kappa\\q):    max dof diff = %.3e\n", diff_sum2);
        if diff_sum2 > 1e-12
            fprintf("Stage 2 sum check FAILED\n"); n_fail = n_fail + 1;
        end
        problem.BoundaryConditions = TransferEssentialBoundaryConditions( ...
            u_, "Dirichlet", radial_map, path_gradient = pg_expr);
    end
    err_u(i) = compute_l2_error(solution.Unknowns.field(2), mesh, u_);
    err_q(i) = compute_l2_error(solution.Unknowns.field(1), mesh, Datum(q_sym));
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

%% -------------------- user gradient reconstruction ------------------- %%

function Gm = darcy_user_gradient(problem, rs, m, n_local, kappa)
% Declares the Darcy constitutive law for the transfer path integral:
%
%   q = -kappa*grad(u)   =>   grad(u).m = -(kappa^{-T} m) . q.
%
% The flux q is local field 1, a vector basis q = sum_i A_i (phi_j q_ij)
% built from a scalar basis phi and component maps A_i (linear_maps), so
% for the basis function phi_j of component block i,
%
%   [G_phi m] = -(kappa^{-T} m)' A_i phi_j.
%
% Inputs (documented path_gradient signature): rs (2 x N x n_b) local
% coordinates of the path points, m (2 x N x n_b) transfer directions.
% Output: Gm (N x n_local x n_b), columns in the local dof layout of
% the local equation block (component blocks of field 1, zero elsewhere).
flux_element = problem.LocalSpaces(1).Basis;
phi = flux_element.ScalarBasis;
d = phi.Dimension;
N = size(rs, 2); n_b = size(rs, 3);

Phi = double(phi.evaluate(reshape(rs, 2, [])));       % (N*n_b x d)
Phi = permute(reshape(Phi, N, n_b, d), [1 3 2]);      % (N x d x n_b)
a = kappa.' \ reshape(m, 2, []);                      % kappa^{-T} m, (2 x N*n_b)

Gm = zeros(N, n_local, n_b);
for i = 1:numel(flux_element.LinearMaps)
    Ai = flux_element.LinearMaps{i};                   % component map (2 x 1)
    ci = -reshape(Ai.' * a, N, 1, n_b);               % -(kappa^{-T} m)' A_i
    Gm(:, (i-1)*d + (1:d), :) = ci .* Phi;
end
end
