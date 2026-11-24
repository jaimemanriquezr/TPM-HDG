% Unfitted Darcy on a level-set domain (BVP-TPM): the kidney-shaped
% region {phi < 0} of KidneyCurve is immersed in a background grid, cut
% elements are discarded, and the Dirichlet datum is transferred from
% the true (implicit) boundary to the polygonal computational boundary
% along closest-point paths (Cockburn & Solano). The curve enters only
% through its level set: no parametrization is ever built. Prints the L2
% error of u with EOC ~ k+1, then plots the mesh and the boundary.
clc; clear;

kappa = eye(2);
k = 2;
kidney = KidneyCurve();

syms x y
u_sym = sin(pi*x) .* sin(pi*y);
q_sym = -kappa * gradient(u_sym, [x y]);
f_sym = divergence(q_sym, [x y]);
u_ = Datum(u_sym);
f_ = Datum(f_sym);

Pk = DubinerBasis(k); VPk = Pk.to_vector(); Mk = LegendreBasis(k);
dx = Measure(preset = "dx");

n_ref = 3;
[err_u, h_list] = deal(zeros(n_ref, 1));
for i = 1:n_ref
    n = 16 * 2^(i-1);
    h_list(i) = 1.6/n;
    background = StructuredSimplicialGrid(n, n, -0.8, 0.8, -0.8, 0.8);
    mesh = immerse_mesh(background, kidney);

    problem = FEProblem(Spaces= [FESpace([VPk, Pk], Condense=true), FESpace(Mk)]);
    [lhs_op, flux_op] = DarcyMatrices(k, mesh, kappa);
    problem.Operator = AffineOperator(lhs_op, flux_op);
    [~, w, ~] = get_test_functions(problem.Space);
    problem.Operator = problem.Operator + assemble(f_ * w * dx, mesh);
    problem.BoundaryConditions = TransferEssentialBoundaryConditions.darcy( ...
        u_, "Dirichlet", kidney, kappa);

    solution = solve(problem, mesh, k);
    err_u(i) = compute_l2_error(solution.Unknowns.field(2), mesh, u_);
    fprintf("n = %3d, h = %.4f, L2 err u = %.4e\n", n, h_list(i), err_u(i));
end
eoc = log(err_u(1:end-1)./err_u(2:end)) ./ log(h_list(1:end-1)./h_list(2:end));
fprintf("EOC(u): %s (expected ~ %d)\n", num2str(eoc.', "%.2f "), k+1);

figure();
plot(mesh); hold on;
plot(kidney, gca(), Color = 'r');
axis equal; title("Unfitted kidney domain: mesh and level-set boundary");
