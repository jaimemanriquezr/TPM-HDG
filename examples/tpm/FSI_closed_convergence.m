% Convergence study of steady FSI on a closed (circular) interface, both
% subdomains meshed fitted to Sigma (fitted_split, closed topology).
% Plots the combined L2 errors and prints EOC; expected k+1 in every
% variable.
clear;
close all;

k = 3;
max_iterations = 6;
h0 = .25;

viscosity = 0.4;
young_modulus = 1.0;
poisson_ratio = 0.2;

uu_fluid = "sin(pi*x)*sin(pi*y)";
p_fluid = "10*x^2";
u_solid = ["sin(pi*x)*cos(pi*y)"; "cos(pi*x)*sin(pi*y)"];

interface = CircularCurve([0.5; 0.5], 0.2, mapping="closest");

[fluid_problem, solid_problem, interface_sources, problem_data] = ...
    FluidStructureCoupling.load_problems(k, interface, 0, 1e-3, ...
        uu_fluid, p_fluid, u_solid, ...
        0, viscosity, 0, young_modulus, poisson_ratio, 0);

[error_sigma, error_u, error_gamma, h_list, n_elements] = deal(zeros(max_iterations, 1));
for i = 1:max_iterations
    disp("Iteration " + i + "/" + max_iterations);
    h1 = h0*(1/2)^(i - 1);
    disp("Generating mesh...")
    [upper_mesh, lower_mesh] = fitted_split(interface, h1, 0);
    disp("Mesh generated.")

    [fluid_solution, solid_solution] = FluidStructureCoupling.solve( ...
        fluid_problem, solid_problem, upper_mesh, lower_mesh, interface, ...
        interface_sources{:});

    error_fluid = fluid_solution.compute_l2_error({problem_data.sigma_F, problem_data.u_F, problem_data.gamma_F}, {problem_data.u_F});
    error_solid = solid_solution.compute_l2_error({problem_data.sigma_S, problem_data.u_S, problem_data.gamma_S}, {problem_data.u_S});

    n_elements(i) = size(upper_mesh.elements, 2) + size(lower_mesh.elements, 2);
    h_list(i) = 1/sqrt(n_elements(i));
    error_sigma(i) = sqrt(error_fluid(1)^2 + error_solid(1)^2);
    error_u(i) = sqrt(error_fluid(2)^2 + error_solid(2)^2);
    error_gamma(i) = sqrt(error_fluid(3)^2 + error_solid(3)^2);
end

r_sigma = compute_eoc(error_sigma, h_list);
r_u = compute_eoc(error_u, h_list);
r_gamma = compute_eoc(error_gamma, h_list);

error_h = @(var) sprintf("$||%s - %s_{h}||_{L^2}$", var, var);
figure();
ax = axes();
loglog(ax, h_list, error_sigma, '*-', 'DisplayName', error_h("\sigma")); hold on;
loglog(ax, h_list, error_u, '*-', 'DisplayName', error_h("u")); hold on;
loglog(ax, h_list, error_gamma, '*-', 'DisplayName', error_h("\gamma")); grid on;

plot_ref = @(p) loglog(h_list, error_u(1) / h_list(1)^(p) / 2 * h_list.^(p), 'k--', 'DisplayName', sprintf("$h^{%i}$", p));
plot_ref(k+1)

legend(ax, 'show', 'Interpreter', 'latex', 'Location', 'best', 'FontSize', 12)
xlabel(ax, "$h$", 'Interpreter', 'latex')
ylabel(ax, "Error", 'Interpreter', 'latex')
title(ax, sprintf('Fitted Stokes/Elasticity --- $\\mathcal{P}_{%i}$', k), 'Interpreter', 'latex', 'FontSize', 16);

e_o_c = [r_sigma, r_u, r_gamma]
