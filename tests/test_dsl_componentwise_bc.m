% ComponentwiseBoundaryConditions on the declarative path must agree
% dof-for-dof with the low-level data.boundary_conditions on the same
% specs. Two problems, both with a mixed tag on x = 1 (outward normal
% (1,0)) fixing one component strongly and the other by its traction:
%   elasticity: u_x dirichlet, (sigma n)_y = sigma_21 neumann
%   fluid:      u_y dirichlet, (sigma n)_x = sigma_11 neumann  (an outflow)
% Polynomial data, so both paths' projections are exact and the tolerance
% can be tight (see test_dsl_elasticity_winkler for why that matters).
clc; clear;
n_fail = 0;
tol = 1e-9;
for k = [1 2]
    Pk = DubinerBasis(k);
    mesh = StructuredSimplicialGrid(8, 8);
    mesh = apply_tag(mesh, "Mixed", @(X) isapprox(X(1, :), 1));

    %% ---- elasticity ----
    E = 1e2; nu = 0.3;
    dS = LinearElasticityProblem.get_manufactured_solution(["x^2*y"; "x*y^2"], E, nu);
    specs = [boundary_condition("Dirichlet", x_data = dS.u{1}, y_data = dS.u{2});
             boundary_condition("Mixed", y = "neumann", x_data = dS.u{1}, y_data = dS.sigma{2, 1})];
    ref = LinearElasticityProblem(k, young_modulus = E, poisson_ratio = nu, source_term = dS.f);
    ref.Data.boundary_conditions = specs;
    sol_ref = solve(ref, mesh, k);

    pb = FEProblem(Order = k, Spaces = [ ...
        FESpace([PEERSBasis(k), Pk.to_vector(), Pk], Condense = true), ...
        FESpace(LegendreBasis(k).to_vector(2))]);
    pb.Mesh = mesh;
    [lhs, flux, rhs] = LinearElasticityMatrices(k, mesh, E, nu);
    pb.Operator = AffineOperator(lhs) + flux;
    if ~isempty(rhs), pb.Operator = pb.Operator + rhs; end
    w = TestFunction(Pk);
    pb.Operator = pb.Operator + IndexedMatrix(IndexedArray(cat(1, ...
        double(Datum(dS.f{1}) * w * Measure(mesh, 'dx')), ...
        double(Datum(dS.f{2}) * w * Measure(mesh, 'dx'))), [-1, 2]));
    pb.BoundaryConditions = ComponentwiseBoundaryConditions(specs);
    sol_dsl = solve(pb, mesh, k);
    d = max(abs(sol_ref.Unknowns.field(2).dof - sol_dsl.Unknowns.field(2).dof), [], 'all');
    e = compute_l2_error(sol_dsl.Unknowns.field(2), mesh, Datum(dS.u));
    fprintf("elasticity k = %d: max|u_dsl - u_lowlevel| = %.3e   e(u) = %.3e\n", k, d, e);
    if d > tol, fprintf(2, "FAILED: elasticity dof mismatch\n"); n_fail = n_fail + 1; end

    %% ---- fluid (incompressible elasticity), pinned pressure ----
    mu = 1.0;
    dF = IncompressibleElasticityProblem.get_manufactured_solution("x^2*y^2", "x^2 + y^2", mu);
    specs = [boundary_condition("Dirichlet", x_data = dF.u{1}, y_data = dF.u{2});
             boundary_condition("Mixed", x = "neumann", x_data = dF.sigma{1, 1}, y_data = dF.u{2})];
    ref = IncompressibleElasticityProblem(k, viscosity = mu, source_term = dF.f, ...
        pressure_mean = "pinned");
    ref.Data.boundary_conditions = specs;
    sol_ref = solve(ref, mesh, k);

    bc = ComponentwiseBoundaryConditions(specs);
    pb = FEProblem(Order = k, Spaces = [ ...
        FESpace([PEERSBasis(k), Pk.to_vector(), Pk, DubinerBasis(0)], Condense = true), ...
        FESpace([LegendreBasis(k).to_vector(2), DubinerBasis(0)])]);
    pb.Mesh = mesh;
    [lhs, flux, rhs] = IncompressibleElasticityMatrices(k, mesh, mu, dF.f, ...
        "pressure_mean", "pinned", natural_faces = bc.natural_faces(mesh));
    pb.Operator = AffineOperator(lhs + rhs) + flux;
    pb.BoundaryConditions = bc;
    sol_dsl = solve(pb, mesh, k);
    d = max(abs(sol_ref.Unknowns.field(2).dof - sol_dsl.Unknowns.field(2).dof), [], 'all');
    e = compute_l2_error(sol_dsl.Unknowns.field(2), mesh, Datum(dF.u));
    fprintf("fluid      k = %d: max|u_dsl - u_lowlevel| = %.3e   e(u) = %.3e\n", k, d, e);
    if d > tol, fprintf(2, "FAILED: fluid dof mismatch\n"); n_fail = n_fail + 1; end
end
fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0, error("test_dsl_componentwise_bc:failed", "%d checks failed", n_fail); end
