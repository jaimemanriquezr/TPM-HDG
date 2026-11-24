% Lambda-robust solid forms: factored condensation and pressure split.
%
% The factored condensation (each element's constant isotropic PEERS stress
% coefficient kept global via FESpace.UncondensedModes) and the pressure
% split (element mean of tr(sigma) kept global, LinearElasticityMatrices
% pressure_split) are exact reformulations for nearly incompressible solids
% and must reproduce the classical discrete displacement to round-off;
% checked on a Dirichlet problem, k = 2, nu = 0.49.
clc; clear;
n_fail = 0;
k = 2; E = 1; nu = 0.49;
Pk = DubinerBasis(k);
mesh = StructuredSimplicialGrid(4, 4);
u_ex = {@(x,y) sin(x).*y; @(x,y) x.^2.*cos(y)};

u = struct();
for variant = ["classical", "factored", "split"]
    if variant == "split"
        spaces = [FESpace([PEERSBasis(k), Pk.to_vector(), Pk, DubinerBasis(0)], Condense = true), ...
            FESpace([LegendreBasis(k).to_vector(2), DubinerBasis(0)])];
    else
        spaces = [FESpace([PEERSBasis(k), Pk.to_vector(), Pk], Condense = true), ...
            FESpace(LegendreBasis(k).to_vector(2))];
    end
    pb = FEProblem(Order = k, Spaces = spaces);
    if variant == "factored"
        pb.Spaces(1).UncondensedModes = {LinearElasticityProblem.isotropic_stress_dofs(k)};
    end
    pb.Mesh = mesh;
    [slhs, sflux, srhs] = LinearElasticityMatrices(k, mesh, E, nu, ...
        "pressure_split", variant == "split");
    pb.Operator = AffineOperator(slhs) + sflux;
    if ~isempty(srhs), pb.Operator = pb.Operator + srhs; end
    w = TestFunction(Pk);
    pb.Operator = pb.Operator + IndexedMatrix(IndexedArray(cat(1, ...
        double(Datum(@(x,y) 1 + 0*x) * w * Measure(mesh, 'dx')), ...
        double(Datum(@(x,y) x.*y) * w * Measure(mesh, 'dx'))), [-1, 2]));
    pb.BoundaryConditions = StrongEssentialBoundaryConditions(Datum(u_ex), ...
        "Dirichlet", projection_order = k);
    sol = solve(pb, mesh, k);
    u.(variant) = sol.Unknowns.field(2).dof;
end

scale = max(abs(u.classical(:)));
d_fac = max(abs(u.factored(:) - u.classical(:))) / scale;
d_spl = max(abs(u.split(:) - u.classical(:))) / scale;
fprintf("relative difference to classical: factored %.2e, split %.2e\n", d_fac, d_spl);
n_fail = n_fail + check(d_fac < 1e-10, "factored condensation reproduces the classical displacement");
n_fail = n_fail + check(d_spl < 1e-10, "pressure split reproduces the classical displacement");

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0, error("test_solid_lambda_forms:failed", "%d checks failed", n_fail); end
