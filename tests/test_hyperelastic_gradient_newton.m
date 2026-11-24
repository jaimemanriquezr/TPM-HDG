%-------------------------------------------------------------------------
%Description: Newton convergence of HyperelasticGradientProblem.solve.
%
%             The manufactured problem of
%             test_hyperelastic_gradient_convergence (u = a (sin(πx) y,
%             x cos(πy)), E = 1, ν = 0.3, Dirichlet data), N = 8, default
%             stabilization, k = 1, 2, from rest (the first step is the
%             linear solve): Saint Venant–Kirchhoff at a = 0.05,
%             neo-Hookean at a = 0.05 and 0.1. Prints |r_n|/|r_0| after
%             each step (newton_residual). Also, the "adaptive"
%             stabilization (τ lagged) at a = 0.05, k = 2, L_c = 0.25.
%             Checks:
%               (1) converged to 1e-10;
%               (2) quadratic convergence: r_{n+1} ≤ 10³ r_n² for
%                   consecutive residuals with r_n < 1e-2 and
%                   r_{n+1} > 1e-13 (default stabilization; the printed
%                   constants are r_{n+1}/r_n²);
%               (3) the stop is the tolerance (stop_reason "tolerance",
%                   stagnated false);
%               (4) with newton_tolerance = 1e-15 (Saint Venant–Kirchhoff,
%                   a = 0.05, k = 1) the stop is the round-off floor:
%                   converged, stagnated, stop_reason "stagnation".
%-------------------------------------------------------------------------
n_failed = 0;
E = 1; nu = 0.3; N = 8;
mesh = StructuredSimplicialGrid(N, N);
cases = {"saint_venant_kirchhoff", 0.05, "mesh_size"; "neo_hookean", 0.05, "mesh_size";
         "neo_hookean", 0.1, "mesh_size"; "saint_venant_kirchhoff", 0.05, "adaptive"};
for c = 1:size(cases, 1)
    [material, a, stabilization] = cases{c, :};
    data = HyperelasticGradientProblem.get_manufactured_solution( ...
        [a + "*sin(pi*x)*y"; a + "*x*cos(pi*y)"], E, nu, 0, 0, 0, material);
    ks = 1:2;
    if stabilization == "adaptive", ks = 2; end
    for k = ks
        pde = HyperelasticGradientProblem(k, material = material, young_modulus = E, poisson_ratio = nu, ...
            source_term = data.f, datum_dirichlet = data.u, stabilization = stabilization, ...
            characteristic_length = 0.25);
        [~, ~, newton] = solve(pde, mesh, k);
        r = newton.residual;
        pairs = find(r(1:end-1) < 1e-2 & r(2:end) > 1e-13);
        constants = r(pairs + 1)./r(pairs).^2;
        fprintf("%s a = %g k = %d %s: |r_n|/|r_0| = %s (r_{n+1}/r_n^2: %s)\n", material, a, k, ...
            stabilization, sprintf("%.1e ", r), sprintf("%.0f ", constants));
        tag = sprintf("%s a = %g k = %d %s: ", material, a, k, stabilization);
        n_failed = n_failed + check(newton.converged && r(end) <= 1e-10, tag + "converged");
        n_failed = n_failed + check(~newton.stagnated && newton.stop_reason == "tolerance", tag + "tolerance stop");
        if stabilization == "mesh_size"
            n_failed = n_failed + check(~isempty(constants) && all(constants <= 1e3), tag + "quadratic");
        end
    end
end
data = HyperelasticGradientProblem.get_manufactured_solution( ...
    ["0.05*sin(pi*x)*y"; "0.05*x*cos(pi*y)"], E, nu, 0, 0, 0, "saint_venant_kirchhoff");
pde = HyperelasticGradientProblem(1, young_modulus = E, poisson_ratio = nu, source_term = data.f, ...
    datum_dirichlet = data.u, newton_tolerance = 1e-15);
[~, ~, newton] = solve(pde, mesh, 1);
fprintf("newton_tolerance = 1e-15: |r_n|/|r_0| = %s(%s)\n", sprintf("%.1e ", newton.residual), newton.stop_reason);
n_failed = n_failed + check(newton.converged && newton.stagnated && newton.stop_reason == "stagnation", ...
    "newton_tolerance = 1e-15: stagnation stop");
fprintf("\n%d check(s) failed.\n", n_failed);
