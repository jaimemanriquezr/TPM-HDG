% test_transfer_bc_elasticity.m
% TPM transfer boundary conditions for PEERS linear elasticity.
%
% Stage 1 (consistency): fitted unit square, identity mapping -> transfer
%   reduces to weak Dirichlet; solution must match dof-for-dof.
% Stage 2 (unfitted): disk R = 0.75 immersed in a background grid of
%   [-1,1]^2, Dirichlet data transferred radially from the circle.
%   Run at nu = 0.3 and at nu = 0.4999 (near-incompressible, the
%   Stokes-stress compliance regime c2 -> 1/(4 mu)).
clc; clear;
n_fail = 0;

k = 2;
E = 1.0;

Pk = DubinerBasis(k);
VPk = Pk.to_vector();
S_k = PEERSBasis(k);
VMk = LegendreBasis(k).to_vector(2);

R = 0.75;
phi = @(X) sqrt(X(1, :, :).^2 + X(2, :, :).^2) - R;
radial_map = @(X) R * X ./ vecnorm(X);

for nu = [0.3, 0.4999]
    lambda = E*nu/((1 + nu)*(1 - 2*nu));
    mu = E/2/(1 + nu);

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

    fprintf("=== nu = %.4f ===\n", nu);

    %% Stage 1: fitted square, identity mapping == weak Dirichlet
    if nu == 0.3
        mesh = StructuredSimplicialGrid(8, 8);
        problem = FEProblem(Spaces= [FESpace([S_k, VPk, Pk], Condense=true), FESpace(VMk)]);
        [lhs_op, flux_op] = LinearElasticityMatrices(k, mesh, E, nu);
problem.Operator = AffineOperator(lhs_op, flux_op);
        w = TestFunction(Pk);
        F1 = double(f1 * w * Measure(mesh, 'dx'));
        F2 = double(f2 * w * Measure(mesh, 'dx'));
        problem.Operator = problem.Operator ...
            + IndexedMatrix(IndexedArray(cat(1, F1, F2), [-1, 2]));

        problem.BoundaryConditions = WeakEssentialBoundaryConditions(u_, "Dirichlet");
        sol_weak = solve(problem, mesh, k);

        problem.BoundaryConditions = TransferEssentialBoundaryConditions.elasticity( ...
            u_, "Dirichlet", @(X) X, E, nu);
        sol_transfer = solve(problem, mesh, k);

        diff_u = max(abs(sol_weak.Unknowns.field(2).dof ...
                       - sol_transfer.Unknowns.field(2).dof), [], 'all');
        fprintf("Stage 1 (identity mapping vs weak BC): max dof diff = %.3e\n", diff_u);
        if diff_u > 1e-9
            fprintf("Stage 1 FAILED\n"); n_fail = n_fail + 1;
        end
    end

    %% Stage 2: unfitted disk, radial transfer
    n_ref = 4;
    [err_u, err_s, h_list] = deal(zeros(n_ref, 1));
    for i = 1:n_ref
        n = 8 * 2^(i-1);
        h_list(i) = 2/n;
        background = StructuredSimplicialGrid(n, n, -1, 1, -1, 1);
        mesh = immerse_mesh(background, phi);

        problem = FEProblem(Spaces= [FESpace([S_k, VPk, Pk], Condense=true), FESpace(VMk)]);
        [lhs_op, flux_op] = LinearElasticityMatrices(k, mesh, E, nu);
problem.Operator = AffineOperator(lhs_op, flux_op);
        w = TestFunction(Pk);
        F1 = double(f1 * w * Measure(mesh, 'dx'));
        F2 = double(f2 * w * Measure(mesh, 'dx'));
        problem.Operator = problem.Operator ...
            + IndexedMatrix(IndexedArray(cat(1, F1, F2), [-1, 2]));
        problem.BoundaryConditions = TransferEssentialBoundaryConditions.elasticity( ...
            u_, "Dirichlet", radial_map, E, nu);

        solution = solve(problem, mesh, k);
        u_h = solution.Unknowns.field(2);
        err_u(i) = compute_l2_error(u_h, mesh, u_);
        fprintf("n = %3d, h = %.4f, L2 err u = %.4e\n", n, h_list(i), err_u(i));
    end
    eoc_u = compute_eoc(err_u, h_list);
    fprintf("EOC(u): %s (expected ~ %d)\n", num2str(eoc_u.', "%.2f "), k+1);

    % rate between the coarsest and finest levels (the mean EOC): the
    % per-level EOCs of the carved disk oscillate around k+1
    rate_u = log(err_u(1)/err_u(end)) / log(h_list(1)/h_list(end));
    fprintf("rate(u) over all levels: %.2f (pass >= %.2f)\n\n", rate_u, k + 1 - 0.25);
    if rate_u < k + 1 - 0.25
        fprintf("Stage 2 rate FAILED (nu = %.4f)\n", nu); n_fail = n_fail + 1;
    end
end

fprintf("%d check(s) failed.\n", n_fail);
