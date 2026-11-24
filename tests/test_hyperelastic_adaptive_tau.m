%-------------------------------------------------------------------------
%Description: The adaptive stabilization of HyperelasticProblem
%             (stabilization = "adaptive", Terrana–Nguyen–Bonet–Peraire
%             eq. 18, τ_F = ρ_max(∂P/∂F)/(2 L_c)).
%
%             Steady manufactured problem u = 0.05 (sin(πx) y, x cos(πy)),
%             E = 1, ν = 0.3, Saint Venant–Kirchhoff, k = 2, N = 8,
%             L_c = 0.5, lagged update.
%             Checks:
%               (1) from rest (no iterate) τ is the linear value
%                   2(μ + λ)/(2 L_c) (the largest eigenvalue of Hooke's law
%                   as a map of the full gradient, plane strain);
%               (2) at the converged iterate τ exceeds the linear value
%                   (the tangent stiffens where the strain grows);
%               (3) Newton converges, and the residual drops by at least
%                   five orders over the last two steps before the
%                   round-off floor (τ is lagged, not differentiated);
%               (4) the u_h error is at most 1.2 times that of the
%                   constant-τ solve with τ = 1.
%-------------------------------------------------------------------------
n_failed = 0;
E = 1; nu = 0.3; k = 2; N = 8; L_c = 0.5;
[lambda, mu] = LinearElasticityProblem.get_lame_parameters(E, nu);
data = HyperelasticProblem.get_manufactured_solution(["0.05*sin(pi*x)*y"; "0.05*x*cos(pi*y)"], E, nu);
mesh = StructuredSimplicialGrid(N, N);
pde = HyperelasticProblem(k, young_modulus = E, poisson_ratio = nu, source_term = data.f, ...
    datum_dirichlet = data.u, stabilization = "adaptive", characteristic_length = L_c);
[solution, ~, newton] = solve(pde, mesh, k);
tau_linear = (mu + lambda)/L_c;
fprintf("tau_max per iterate: %s (linear value %.4f)\n", sprintf("%.4f ", newton.tau_max), tau_linear);
fprintf("relative residuals: %s\n", sprintf("%.2e ", newton.residual));
n_failed = n_failed + check(abs(newton.tau_max(1) - tau_linear) <= 1e-12*tau_linear, "linear tau from rest");
n_failed = n_failed + check(newton.tau_max(end) > tau_linear, "tau grows with the strain");
r = newton.residual;
n_failed = n_failed + check(newton.converged && r(end) <= 1e-10, "Newton converged");
n_failed = n_failed + check(numel(r) >= 3 && r(end - 1)/r(end - 3) <= 1e-5, ...
    "fast Newton convergence with the lagged tau");

constant = HyperelasticProblem(k, young_modulus = E, poisson_ratio = nu, source_term = data.f, ...
    datum_dirichlet = data.u);
e_adaptive = solution.compute_l2_error({data.sigma, data.u, data.gamma}, {data.u});
e_constant = solve(constant, mesh, k).compute_l2_error({data.sigma, data.u, data.gamma}, {data.u});
fprintf("e(u_h): adaptive %.3e, constant tau = 1 %.3e\n", e_adaptive(2), e_constant(2));
n_failed = n_failed + check(e_adaptive(2) <= 1.2*e_constant(2), ...
    "u_h error no larger than with constant tau");
fprintf("\n%d check(s) failed.\n", n_failed);
