% Convergence study of steady FSI on a smooth (sinusoidal) interface,
% both subdomains meshed fitted to Sigma with no gap. Plots the errors and
% prints EOC; expected k+1 in every variable.
clear;
close all;

k = 3;
max_iterations = 4;
h0 = .25 / 2;

viscosity = 0.4;
young_modulus = 1.0;
poisson_ratio = 0.2;

uu_fluid = "1/pi*sin(pi*x)*sin(pi*y)";
p_fluid = "10*x^2";
u_solid = ["-sin(pi*x)*cos(pi*y)"; "cos(pi*x)*sin(pi*y)"];

interface = SmoothCurveHorizontal(mapping="closest");

[fluid_problem, solid_problem, interface_sources, problem_data] = ...
    FluidStructureCoupling.load_problems(k, interface, 0, 1e-3, ...
        uu_fluid, p_fluid, u_solid, ...
        0, viscosity, 0, young_modulus, poisson_ratio, 0);

[error_sigma, error_u, error_gamma, h_list, n_elements] = deal(zeros(max_iterations, 1));
for i = 1:max_iterations
    disp("Iteration " + i + "/" + max_iterations);
    h1 = h0*(3/4)^(i + 1);
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

figure();
loglog(h_list, error_sigma, '*-', 'DisplayName', '\sigma'); hold on
loglog(h_list, error_u, '*-', 'DisplayName', 'u')
loglog(h_list, error_gamma, '*-', 'DisplayName', '\gamma')
grid on
legend show

e_o_c = [r_sigma, r_u, r_gamma]
mean(e_o_c)
