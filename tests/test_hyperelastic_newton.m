%-------------------------------------------------------------------------
%Description: Newton convergence of HyperelasticProblem.solve.
%
%             The steady manufactured problem of test_hyperelastic_convergence
%             (u = a (sin(πx) y, x cos(πy)), E = 1, ν = 0.3, Dirichlet
%             data) at a = 0.05 (|∇u| up to 0.16), N = 8, τ = 4/h, both
%             materials and k = 1, 2, from the zero iterate (the first step
%             is the linear solve). Prints the relative residual
%             |r_n| / |r_0| after each step (newton_residual).
%             Checks:
%               (1) converged to the tolerance 1e-10;
%               (2) quadratic convergence: for consecutive residuals with
%                   1e-13 < r_{n+1} and r_n < 1e-3, log r_{n+1} / log r_n ≥ 1.5
%                   (2 for exact quadratic convergence);
%               (3) the stop is the tolerance (stop_reason "tolerance",
%                   stagnated false);
%               (4) with newton_tolerance = 1e-15 (Saint Venant–Kirchhoff,
%                   k = 1) the stop is the round-off floor: converged,
%                   stagnated, stop_reason "stagnation".
%-------------------------------------------------------------------------
n_failed = 0;
E = 1; nu = 0.3; N = 8;
u_exact = ["0.05*sin(pi*x)*y"; "0.05*x*cos(pi*y)"];
mesh = StructuredSimplicialGrid(N, N).add_stabilization(4*N);
for material = ["saint_venant_kirchhoff", "neo_hookean"]
    data = HyperelasticProblem.get_manufactured_solution(u_exact, E, nu, 0, 0, 0, material);
    for k = 1:2
        pde = HyperelasticProblem(k, material = material, young_modulus = E, poisson_ratio = nu, ...
            source_term = data.f, datum_dirichlet = data.u);
        [solution, ~, newton] = solve(pde, mesh, k);
        r = newton.residual;
        fprintf("%s, k = %d: |r_n|/|r_0| = %s\n", material, k, sprintf("%.2e  ", r));
        tag = sprintf("%s k = %d: ", material, k);
        n_failed = n_failed + check(newton.converged && r(end) <= 1e-10, tag + "converged");
        pairs = find(r(1:end-1) < 1e-3 & r(2:end) > 1e-13);
        orders = log(r(pairs + 1))./log(r(pairs));
        fprintf("    observed orders log r_{n+1}/log r_n: %s\n", sprintf("%.2f  ", orders));
        n_failed = n_failed + check(~isempty(orders) && all(orders >= 1.5), tag + "quadratic convergence");
        n_failed = n_failed + check(~newton.stagnated && newton.stop_reason == "tolerance", tag + "tolerance stop");
    end
end
data = HyperelasticProblem.get_manufactured_solution(u_exact, E, nu, 0, 0, 0, "saint_venant_kirchhoff");
pde = HyperelasticProblem(1, young_modulus = E, poisson_ratio = nu, source_term = data.f, ...
    datum_dirichlet = data.u, newton_tolerance = 1e-15);
[~, ~, newton] = solve(pde, mesh, 1);
fprintf("newton_tolerance = 1e-15: |r_n|/|r_0| = %s(%s)\n", sprintf("%.2e  ", newton.residual), newton.stop_reason);
n_failed = n_failed + check(newton.converged && newton.stagnated && newton.stop_reason == "stagnation", ...
    "newton_tolerance = 1e-15: stagnation stop");
fprintf("\n%d check(s) failed.\n", n_failed);
