% test_transfer_bc_elasticity_declarative.m
% Declarative path_gradient for PEERS linear elasticity:
%
%   pg = c1*sigma - c2*trace(sigma)*eye(2) + p*[0 1; -1 0]
%
% lowered through expression_gradient's mesh-dependent branch
% (PEERSBasis.evaluate per-page mode), verified dof-for-dof against the
% .elasticity factory -- an independent hand-computed implementation
% (elasticity_gradient) -- then unfitted-disk EOC at nu = 0.3 and 0.4999.
%
% Stage 1 (fitted square, identity mapping) has m = 0, so it checks the
% plumbing only; stage 2's first refinement is the real lowering check.
clc; clear;

k = 2;
E = 1.0;

Pk = DubinerBasis(k);
VPk = Pk.to_vector();
S_k = PEERSBasis(k);
VMk = LegendreBasis(k).to_vector(2);

R = 0.75;
phi = @(X) sqrt(X(1, :, :).^2 + X(2, :, :).^2) - R;
radial_map = @(X) R * X ./ vecnorm(X);

n_fail = 0;
for nu = [0.3, 0.4999]
    lambda = E*nu/((1 + nu)*(1 - 2*nu));
    mu = E/2/(1 + nu);
    c1 = 1/(2*mu);
    c2 = lambda/(4*mu*(lambda + mu));

    syms x y
    u_sym = [sin(pi*x) .* sin(pi*y); cos(pi*x) .* cos(pi*y)];
    gu = [diff(u_sym(1), x), diff(u_sym(1), y);
          diff(u_sym(2), x), diff(u_sym(2), y)];
    eps_u = (gu + gu.')/2;
    sig = lambda*trace(eps_u)*eye(2) + 2*mu*eps_u;
    f_sym = -[diff(sig(1,1), x) + diff(sig(1,2), y);
              diff(sig(2,1), x) + diff(sig(2,2), y)];
    u_ = Datum(u_sym);
    f1 = Datum(f_sym(1));
    f2 = Datum(f_sym(2));

    fprintf("=== nu = %.4f (c1 = %.4f, c2 = %.4f) ===\n", nu, c1, c2);

    %% Stage 1: fitted square, identity mapping, declarative vs factory
    if nu == 0.3
        mesh = StructuredSimplicialGrid(8, 8);
        problem = build_problem(mesh, S_k, VPk, Pk, VMk, k, E, nu, f1, f2);
        [sigma, ~, p] = get_trial_functions(problem.Space);
        pg = c1*sigma - c2*trace(sigma)*eye(2) + p*[0 1; -1 0];

        problem.BoundaryConditions = TransferEssentialBoundaryConditions.elasticity( ...
            u_, "Dirichlet", @(X) X, E, nu);
        sol_factory = solve(problem, mesh, k);

        problem.BoundaryConditions = TransferEssentialBoundaryConditions( ...
            u_, "Dirichlet", @(X) X, path_gradient = pg);
        sol_expr = solve(problem, mesh, k);

        diff_u = max(abs(sol_factory.Unknowns.field(2).dof ...
                       - sol_expr.Unknowns.field(2).dof), [], 'all');
        fprintf("Stage 1 (declarative vs factory, m = 0): max dof diff = %.3e\n", diff_u);
        if diff_u > 1e-12
            fprintf("Stage 1 FAILED\n"); n_fail = n_fail + 1;
        end
    end

    %% Stage 2: unfitted disk, radial transfer, declarative BC
    n_ref = 4;
    [err_u, h_list] = deal(zeros(n_ref, 1));
    for i = 1:n_ref
        n = 8 * 2^(i-1);
        h_list(i) = 2/n;
        background = StructuredSimplicialGrid(n, n, -1, 1, -1, 1);
        mesh = immerse_mesh(background, phi);

        problem = build_problem(mesh, S_k, VPk, Pk, VMk, k, E, nu, f1, f2);
        [sigma, ~, p] = get_trial_functions(problem.Space);
        pg = c1*sigma - c2*trace(sigma)*eye(2) + p*[0 1; -1 0];
        problem.BoundaryConditions = TransferEssentialBoundaryConditions( ...
            u_, "Dirichlet", radial_map, path_gradient = pg);

        solution = solve(problem, mesh, k);
        u_h = solution.Unknowns.field(2);
        err_u(i) = compute_l2_error(u_h, mesh, u_);
        fprintf("n = %3d, h = %.4f, L2 err u = %.4e\n", n, h_list(i), err_u(i));

        if i == 1
            % the real lowering check: nonzero m, factory as the oracle
            problem.BoundaryConditions = TransferEssentialBoundaryConditions.elasticity( ...
                u_, "Dirichlet", radial_map, E, nu);
            sol_factory = solve(problem, mesh, k);
            diff_u = max(abs(sol_factory.Unknowns.field(2).dof - u_h.dof), [], 'all');
            fprintf("Stage 2 lvl 1 (declarative vs factory): max dof diff = %.3e\n", diff_u);
            if diff_u > 1e-12
                fprintf("Stage 2 dof comparison FAILED\n"); n_fail = n_fail + 1;
            end
        end
    end
    eoc_u = compute_eoc(err_u, h_list);
    fprintf("EOC(u): %s (expected ~ %d)\n\n", num2str(eoc_u.', "%.2f "), k+1);
    % the per-level rate wobbles on the unfitted disk (same trajectory as
    % the factory test); the mean over the ladder is the stable statistic
    if mean(eoc_u) < k + 0.5
        fprintf("EOC check FAILED\n"); n_fail = n_fail + 1;
    end
end
fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_transfer_bc_elasticity_declarative:failed", "%d checks failed", n_fail);
end

function problem = build_problem(mesh, S_k, VPk, Pk, VMk, k, E, nu, f1, f2)
problem = FEProblem(Spaces= [FESpace([S_k, VPk, Pk], Condense=true), FESpace(VMk)]);
[lhs_op, flux_op] = LinearElasticityMatrices(k, mesh, E, nu);
problem.Operator = AffineOperator(lhs_op, flux_op);
w = TestFunction(Pk);
F1 = double(f1 * w * Measure(mesh, 'dx'));
F2 = double(f2 * w * Measure(mesh, 'dx'));
problem.Operator = problem.Operator ...
    + IndexedMatrix(IndexedArray(cat(1, F1, F2), [-1, 2]));
end
