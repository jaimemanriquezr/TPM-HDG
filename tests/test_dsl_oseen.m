% test_dsl_oseen.m
% Oseen / Navier-Stokes on the declarative path (StokesMatrices /
% IncompressibleElasticityMatrices with a convective velocity, as
% IndexedMatrix blocks) against the low-level classes, dof-for-dof:
%   (1) fitted square, prescribed beta, weak Dirichlet == strong Dirichlet
%       of OseenProblem / NavierStokesProblem (the datum is polynomial so
%       both projections coincide);
%   (2) immersed disk, radial TPM transfer: the declarative
%       TransferEssentialBoundaryConditions (velocity_gradient /
%       incompressible_stress factories) against the low-level
%       transfer_boundary_rows -- the two TPM implementations check each
%       other, with convection on;
%   (3) declarative Picard Navier-Stokes (operator_of) on the disk against
%       the low-level Picard.
clc; clear;
n_fail = 0;
nu = 0.1;
beta = ["1 + y", "x"];
R = 0.75;
phi = @(X) sqrt(X(1, :, :).^2 + X(2, :, :).^2) - R;
radial = @(X) R * X ./ vecnorm(X);

Pk_of = @(k) DubinerBasis(k);
spaces_L = @(k) [FESpace([Pk_of(k).to_tensor(), Pk_of(k).to_vector(), Pk_of(k), DubinerBasis(0)], Condense=true), ...
                 FESpace([LegendreBasis(k).to_vector(2), DubinerBasis(0), ConstantBasis()])];
spaces_S = @(k) [FESpace([PEERSBasis(k), Pk_of(k).to_vector(), Pk_of(k), DubinerBasis(0)], Condense=true), ...
                 FESpace([LegendreBasis(k).to_vector(2), DubinerBasis(0), ConstantBasis()])];

for stage = 1:2
    if stage == 1
        U = "x^2 + x*y + y^2";  P = "x*y";  meshes = {StructuredSimplicialGrid(8, 8)};  label = "fitted square";
    else
        U = "sin(pi*x)*sin(pi*y)";  P = "cos(pi*x)*sin(pi*y)";
        meshes = {immerse_mesh(StructuredSimplicialGrid(16, 16, -1, 1, -1, 1), phi)};  label = "immersed disk (transfer)";
    end
    dL = OseenProblem.get_manufactured_solution(U, P, nu, beta);
    dS = NavierStokesProblem.get_manufactured_solution(U, P, nu, beta);
    mesh = meshes{1};
    for k = [1 2]
        %% velocity-gradient form
        pL = OseenProblem(k, viscosity = nu, source_term = dL.f, datum_dirichlet = dL.u, convective_velocity = dL.beta);
        if stage == 2, pL.Data.transfer_mapping = radial; end
        sol2 = solve(pL, mesh, Order = k);
        problem = FEProblem(Spaces = spaces_L(k));
        problem.Mesh = mesh;
        [lhs, flux, rhs] = StokesMatrices(k, mesh, nu, dL.f, dL.beta);
        problem.Operator = AffineOperator(lhs + rhs) + flux;
        if stage == 1
            problem.BoundaryConditions = WeakEssentialBoundaryConditions(Datum(dL.u), "Dirichlet");
        else
            problem.BoundaryConditions = TransferEssentialBoundaryConditions.velocity_gradient(Datum(dL.u), "Dirichlet", radial);
        end
        sol1 = solve(problem, mesh, k);
        dfL = max(abs(sol1.Unknowns.field(2).dof - sol2.Unknowns.field(2).dof), [], 'all');
        dgL = max(abs(sol1.Unknowns.field(1).dof - sol2.Unknowns.field(1).dof), [], 'all');

        %% stress form
        pS = oseen_stress(dS.beta, k, viscosity = nu, source_term = dS.f, datum_dirichlet = dS.u);
        if stage == 2, pS.Data.transfer_mapping = radial; end
        sol2 = solve(pS, mesh, Order = k);
        problem = FEProblem(Spaces = spaces_S(k));
        problem.Mesh = mesh;
        [lhs, flux, rhs] = IncompressibleElasticityMatrices(k, mesh, nu, dS.f, convective_velocity = dS.beta);
        problem.Operator = AffineOperator(lhs + rhs) + flux;
        if stage == 1
            problem.BoundaryConditions = WeakEssentialBoundaryConditions(Datum(dS.u), "Dirichlet");
        else
            problem.BoundaryConditions = TransferEssentialBoundaryConditions.incompressible_stress(Datum(dS.u), "Dirichlet", radial, nu);
        end
        sol1 = solve(problem, mesh, k);
        dfS = max(abs(sol1.Unknowns.field(2).dof - sol2.Unknowns.field(2).dof), [], 'all');
        dgS = max(abs(sol1.Unknowns.field(1).dof - sol2.Unknowns.field(1).dof), [], 'all');

        scale = max(abs(sol2.Unknowns.field(2).dof), [], 'all');
        fprintf("Stage %d (%s) k = %d: declarative vs low-level  L-form |du| %.1e |dL| %.1e   S-form |du| %.1e |dsigma| %.1e  (|u| ~ %.1e)\n", ...
            stage, label, k, dfL, dgL, dfS, dgS, scale);
        if max([dfL, dgL, dfS, dgS]) > 1e-9 * max(1, scale)
            fprintf("Stage %d k = %d FAILED\n", stage, k); n_fail = n_fail + 1;
        end
    end
end

%% Stage 3: declarative Picard on the disk (velocity-gradient form)
k = 2;
mesh = immerse_mesh(StructuredSimplicialGrid(16, 16, -1, 1, -1, 1), phi);
dL = OseenProblem.get_manufactured_solution("sin(pi*x)*sin(pi*y)", "cos(pi*x)*sin(pi*y)", nu, []);
pL = OseenProblem(k, viscosity = nu, source_term = dL.f, datum_dirichlet = dL.u, transfer_mapping = radial);
[s2, h2] = solve_navier_stokes(pL, mesh, k, tolerance = 1e-10);
problem = FEProblem(Spaces = spaces_L(k));
problem.Mesh = mesh;
problem.Data.viscosity = nu;
problem.BoundaryConditions = TransferEssentialBoundaryConditions.velocity_gradient(Datum(dL.u), "Dirichlet", radial);
operator_of = @(b) oseen_operator(k, mesh, nu, dL.f, b);
[s1, h1] = solve_navier_stokes(problem, mesh, k, tolerance = 1e-10, operator_of = operator_of);
du = max(abs(s1.Unknowns.field(2).dof - s2.Unknowns.field(2).dof), [], 'all');
fprintf("Stage 3 Picard: declarative %d it, low-level %d it, |du| = %.1e\n", numel(h1), numel(h2), du);
if du > 1e-8 || numel(h1) ~= numel(h2)
    fprintf("Stage 3 FAILED\n"); n_fail = n_fail + 1;
end

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_dsl_oseen:failed", "%d checks failed", n_fail);
end

function op = oseen_operator(k, mesh, nu, f, beta)
[lhs, flux, rhs] = StokesMatrices(k, mesh, nu, f, beta);
op = AffineOperator(lhs + rhs) + flux;
end
