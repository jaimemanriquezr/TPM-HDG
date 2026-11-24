%-------------------------------------------------------------------------
%Description: NavierStokesProblem.solve is always a Picard iteration
%             (Stokes, or the problem's convective velocity as the starting
%             beta, then beta = u* of the previous iterate). It reproduces a
%             hand-written loop over the parent's single solves; a plain
%             FEProblem with the class's factory blocks, iterated by
%             solve_navier_stokes (operator_of), agrees with it; and a
%             restart from the converged u* reconverges at once. Convective
%             stabilization is off throughout (nu = 1, so Picard converges
%             without it).
%-------------------------------------------------------------------------
n_failed = 0;
nu = 1;  k = 2;  tol = 1e-10;
U = "sin(pi*x)^2*sin(pi*y)^2";  P = "cos(pi*x)*sin(pi*y)";
mesh = StructuredSimplicialGrid(8, 8);
d = NavierStokesProblem.get_manufactured_solution(U, P, nu, []);

%% reference: a hand-written Picard loop over single solves --------------
ref = IncompressibleElasticityProblem(k, viscosity = nu, source_term = d.f, ...
    datum_dirichlet = d.u, convective_stabilization = false);
beta = {};  u_old = [];
for n = 1:50
    ref.Data.convective_velocity = beta;
    s_ref = solve(ref, mesh, Order = k);
    u = s_ref.Unknowns.field(2).dof;
    beta = divergence_free_velocity(s_ref);
    if ~isempty(u_old) && norm(u - u_old, "fro")/norm(u_old, "fro") < tol
        break
    end
    u_old = u;
end
u_ref = s_ref.Unknowns.field(2).dof;
scale = max(abs(u_ref), [], "all");

%% low-level path ---------------------------------------------------------
p = NavierStokesProblem(k, viscosity = nu, source_term = d.f, datum_dirichlet = d.u, ...
    convective_stabilization = false);
[s1, ~, ~, u_star] = solve(p, mesh, k, tolerance = tol);
du = max(abs(s1.Unknowns.field(2).dof - u_ref), [], "all");
n_failed = n_failed + check(sprintf("low-level vs hand-written loop: |du| = %.1e", du), ...
    du <= 1e-12*scale);

%% a plain FEProblem with the class's factory blocks (solve_navier_stokes, operator_of)
Pk = DubinerBasis(k);
pd = FEProblem(Spaces = [FESpace([PEERSBasis(k), Pk.to_vector(), Pk, DubinerBasis(0)], Condense = true), ...
                         FESpace([LegendreBasis(k).to_vector(2), DubinerBasis(0), ConstantBasis()])]);
pd.Mesh = mesh;
pd.Data.viscosity = nu;   % read by divergence_free_velocity
pd.BoundaryConditions = WeakEssentialBoundaryConditions(Datum(d.u), "Dirichlet");
s2 = solve_navier_stokes(pd, mesh, k, tolerance = tol, operator_of = @(b) factory_operator(k, mesh, nu, d.f, b));
du = max(abs(s2.Unknowns.field(2).dof - u_ref), [], "all");
n_failed = n_failed + check(sprintf("plain FEProblem (factory operator) vs low-level: |du| = %.1e", du), ...
    du <= 1e-9*scale);

%% a convective velocity on the problem is the starting beta of the Picard loop
p.Data.convective_velocity = u_star;
[s3, ~, h3] = solve(p, mesh, k, tolerance = tol);
du = max(abs(s3.Unknowns.field(2).dof - u_ref), [], "all");
n_failed = n_failed + check(sprintf("restart from the converged u*: %d Oseen step(s), |du| = %.1e", ...
    numel(h3), du), numel(h3) <= 2 && du <= 1e-9*scale);

fprintf("\n%d check(s) failed.\n", n_failed);

%=========================================================================%
function failed = check(name, condition)
failed = ~condition;
if failed
    fprintf(2, "  FAIL  %s\n", name);
else
    fprintf("  ok    %s\n", name);
end
end

function op = factory_operator(k, mesh, nu, f, beta)
opts = {};
if ~isempty(beta), opts = {"convective_velocity", beta, "convective_stabilization", false}; end
[lhs, flux, rhs] = IncompressibleElasticityMatrices(k, mesh, nu, f, opts{:});
op = AffineOperator(lhs + rhs) + flux;
end
