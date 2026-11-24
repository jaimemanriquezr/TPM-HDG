% test_transfer_bc_levelset.m
% A LevelSetCurve as the unfitted BVP boundary (BVP-TPM): the disk of
% test_transfer_bc_darcy stage 2, once with the hand-built phi/radial-map
% handles (the established oracle) and once passing the LevelSetCurve
% itself to immerse_mesh and the transfer BC. The curve's Newton
% projection must reproduce the radial map, so the two solves agree
% dof-for-dof up to the projection tolerance.
clc; clear;

kappa = eye(2);
k = 2;
R = 0.75;
n = 16;
n_fail = 0;

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

background = StructuredSimplicialGrid(n, n, -1, 1, -1, 1);

% oracle: hand-built handles
phi = @(X) sqrt(X(1, :, :).^2 + X(2, :, :).^2) - R;
radial_map = @(X) R * X ./ vecnorm(X);

% the same disk as a LevelSetCurve (analytic gradient)
circle = LevelSetCurve(@(X) X(1, :).^2 + X(2, :).^2 - R^2, ...
    gradient = @(X) 2*X, box = [-1 1 -1 1]);

solutions = cell(1, 2);
geometries = {{phi, radial_map}, {circle, circle}};
for i = 1:2
    [carve, mapping] = geometries{i}{:};
    mesh = immerse_mesh(background, carve);
    problem = FEProblem(Spaces= [FESpace([VPk, Pk], Condense=true), FESpace(Mk)]);
    [lhs_op, flux_op] = DarcyMatrices(k, mesh, kappa);
    problem.Operator = AffineOperator(lhs_op, flux_op);
    [~, w, ~] = get_test_functions(problem.Space);
    problem.Operator = problem.Operator + assemble(f_ * w * dx, mesh);
    problem.BoundaryConditions = TransferEssentialBoundaryConditions.darcy( ...
        u_, "Dirichlet", mapping, kappa);
    solutions{i} = solve(problem, mesh, k);
end

diff_u = max(abs(solutions{1}.Unknowns.field(2).dof - solutions{2}.Unknowns.field(2).dof), [], 'all');
diff_q = max(abs(solutions{1}.Unknowns.field(1).dof - solutions{2}.Unknowns.field(1).dof), [], 'all');
fprintf("LevelSetCurve vs handle oracle: max dof diff u = %.3e, q = %.3e\n", diff_u, diff_q);
if diff_u > 1e-8 || diff_q > 1e-8
    fprintf("FAIL level-set transfer disagrees with handle oracle\n");
    n_fail = n_fail + 1;
end

err_u = compute_l2_error(solutions{2}.Unknowns.field(2), ...
    immerse_mesh(background, circle), u_);
fprintf("L2 err u (level-set run) = %.4e\n", err_u);
if err_u > 1e-2
    fprintf("FAIL level-set solve error out of range\n");
    n_fail = n_fail + 1;
end

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_transfer_bc_levelset:failed", "%d checks failed", n_fail);
end
