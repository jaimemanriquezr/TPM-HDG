% test_biot_componentwise_bc.m
% Componentwise boundary conditions on the three-component Biot trace
% (eta_x, eta_y, p): boundary_condition(..., p = ...) and
% ComponentwiseBoundaryConditions.
%
%   (a) a spec with x, y and p "dirichlet" on every face equals the class
%       default (strong Dirichlet data {eta_x; eta_y; p} on "Dirichlet")
%       dof for dof;
%   (b) a two-component spec leaves p_hat natural with zero data: it
%       equals the three-component spec with p = "neumann" and zero data
%       exactly, and differs from (a);
%   (c) a spec with more components than the trace errors.
% Unit square, 4 x 4, k = 2, one Euler step from rest.
clc; clear;
n_fail = 0;

k = 2;
mesh = StructuredSimplicialGrid(4, 4);       % every boundary face "Dirichlet"
gx = @(x, y) x.^2 + 0*y;  gy = @(x, y) x.*y;  gp = @(x, y) x + 2*y;
problem = BiotProblem(k, young_modulus = 1, poisson_ratio = 0.3, biot_alpha = 1, ...
    storage = 1, diffusion = eye(2), time_step = 0.5, ...
    source_term = {@(x, y) 1 + 0*x; @(x, y) x - y}, source_mass = @(x, y) 1 + x.*y, ...
    datum_dirichlet = {gx; gy; gp});
trace = @(s) s.Unknowns.hybrid(1).dof;

%% (a) three dirichlet components = the strong default
reference = solve(problem, mesh, k);
problem.Data.boundary_conditions = boundary_condition("Dirichlet", ...
    x_data = gx, y_data = gy, p = "dirichlet", p_data = gp);
spec_3 = solve(problem, mesh, k);
err = max(abs(trace(spec_3) - trace(reference)), [], "all");
fprintf("(a) max trace difference to the strong default: %.2e\n", err);
n_fail = n_fail + check(err < 1e-10, "(a) p = ""dirichlet"" equals the strong default");

%% (b) two components: p_hat natural, zero flux
problem.Data.boundary_conditions = boundary_condition("Dirichlet", x_data = gx, y_data = gy);
spec_2 = solve(problem, mesh, k);
problem.Data.boundary_conditions = boundary_condition("Dirichlet", ...
    x_data = gx, y_data = gy, p = "neumann");
spec_n = solve(problem, mesh, k);
n_fail = n_fail + check(isequal(trace(spec_2), trace(spec_n)), ...
    "(b) no p equals p = ""neumann"" with zero data");
gap = max(abs(trace(spec_2) - trace(reference)), [], "all");
fprintf("(b) max trace difference to (a): %.2e\n", gap);
n_fail = n_fail + check(gap > 1e-6, "(b) p_hat is not prescribed by a two-component spec");

%% (c) too many components
bad = boundary_condition("Dirichlet", p = "dirichlet");
bad.type = [bad.type; "dirichlet"];
bad.data = [bad.data; {gp}];
problem.Data.boundary_conditions = bad;
try
    solve(problem, mesh, k);
    raised = "";
catch exception
    raised = string(exception.identifier);
end
n_fail = n_fail + check(raised == "ComponentwiseBoundaryConditions:components", ...
    "(c) four components on a three-component trace error");

fprintf("%d check(s) failed.\n", n_fail);

function failed = check(ok, label)
if ok
    fprintf("  pass  %s\n", label);
else
    fprintf("  FAIL  %s\n", label);
end
failed = double(~ok);
end
