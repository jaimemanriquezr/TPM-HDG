% Winkler foundation on the declarative elasticity builder.
% LinearElasticityMatrices(..., hooke_modulus = b0) must agree dof-for-dof
% with the low-level LinearElasticityProblem(hooke_modulus = b0). The
% manufactured source carries the b0*u term, so the exact solution is the
% same on both paths.
clc; clear;
n_fail = 0;
E = 1e2; nu = 0.3; b0 = 4e2;
for k = [1 2]
    % data of degree k+1: the Winkler source b0*u then integrates exactly
    % against degree-k test functions on both paths
    if k == 1, u_exact = ["x + 2*y"; "3*x - y"]; else, u_exact = ["x^2*y"; "x*y^2"]; end
    data = LinearElasticityProblem.get_manufactured_solution(u_exact, E, nu, b0);
    Pk = DubinerBasis(k);
    mesh = StructuredSimplicialGrid(8, 8);
    % --- low-level reference ---
    ref = LinearElasticityProblem(k, young_modulus = E, poisson_ratio = nu, ...
        hooke_modulus = b0, source_term = data.f, datum_dirichlet = data.u);
    sol_ref = solve(ref, mesh, k);
    % --- declarative ---
    pb = FEProblem(Order = k, Spaces = [ ...
        FESpace([PEERSBasis(k), Pk.to_vector(), Pk], Condense = true), ...
        FESpace(LegendreBasis(k).to_vector(2))]);
    pb.Mesh = mesh;
    [lhs, flux, rhs] = LinearElasticityMatrices(k, mesh, E, nu, hooke_modulus = b0);
    pb.Operator = AffineOperator(lhs) + flux;
    if ~isempty(rhs), pb.Operator = pb.Operator + rhs; end
    w = TestFunction(Pk);
    pb.Operator = pb.Operator + IndexedMatrix(IndexedArray(cat(1, ...
        double(Datum(data.f{1}) * w * Measure(mesh, 'dx')), ...
        double(Datum(data.f{2}) * w * Measure(mesh, 'dx'))), [-1, 2]));
    % the same skeleton projection the low-level path uses, so the only
    % difference left between the two paths is the Winkler term itself
    pb.BoundaryConditions = ComponentwiseBoundaryConditions( ...
        boundary_condition("Dirichlet", x_data = data.u{1}, y_data = data.u{2}));
    sol_dsl = solve(pb, mesh, k);
    d = max(abs(sol_ref.Unknowns.field(2).dof - sol_dsl.Unknowns.field(2).dof), [], 'all');
    e = compute_l2_error(sol_dsl.Unknowns.field(2), mesh, Datum(data.u));
    fprintf("k = %d: max|u_dsl - u_lowlevel| = %.3e   e(u) = %.3e\n", k, d, e);
    if d > 1e-10, fprintf(2, "FAILED: dof mismatch\n"); n_fail = n_fail + 1; end
    if e > 5e-2,  fprintf(2, "FAILED: not converged\n"); n_fail = n_fail + 1; end
end
fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0, error("test_dsl_elasticity_winkler:failed", "%d checks failed", n_fail); end
