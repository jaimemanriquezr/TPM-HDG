% Convergence study of steady FSI on a flat interface with a gap closing
% as h^2 (fitted meshes retreating from y0; TPM bridges the band).
% Prints EOC; expected k+1 in every variable.
clc; clear;
close all;

k = 2;
max_iterations = 5;
h0 = .3;

y0 = .5;
gap = @(h) h^2;

viscosity = 1.00;
young_modulus = 1.0;
poisson_ratio = 0.33;

uu_fluid = "sin(pi*x)*sin(pi*y)";
p_fluid = "0";
u_solid = ["cos(pi*x)"; "sin(pi*y)"];

interface = FlatCurveHorizontal(y0);

[fluid_problem, solid_problem, interface_sources, problem_data] = ...
    FluidStructureCoupling.load_problems(k, interface, 0, 1e-3, ...
        uu_fluid, p_fluid, u_solid, ...
        0, viscosity, 0, young_modulus, poisson_ratio, 0);

[error_sigma, error_u, error_gamma, error_hat, h_list, n_elements] = deal(zeros(max_iterations, 1));
for i = 1:max_iterations
    disp("Iteration " + i + "/" + max_iterations);
    h1 = h0*(1/2)^(i - 1);
    disp("Generating mesh...")
    [upper_mesh, lower_mesh] = fitted_split(interface, h1, gap(h1));
    disp("Mesh generated.")

    [fluid_solution, solid_solution] = FluidStructureCoupling.solve( ...
        fluid_problem, solid_problem, upper_mesh, lower_mesh, interface, ...
        interface_sources{:});

    [error_fluid, error_hat_fluid] = fluid_solution.compute_l2_error({problem_data.sigma_F, problem_data.u_F, problem_data.gamma_F}, {problem_data.u_F});
    [error_solid, error_hat_solid] = solid_solution.compute_l2_error({problem_data.sigma_S, problem_data.u_S, problem_data.gamma_S}, {problem_data.u_S});

    n_elements(i) = size(upper_mesh.elements, 2) + size(lower_mesh.elements, 2);
    h_list(i) = 1/sqrt(n_elements(i));
    error_sigma(i) = sqrt(error_fluid(1)^2 + error_solid(1)^2);
    error_u(i) = sqrt(error_fluid(2)^2 + error_solid(2)^2);
    error_gamma(i) = sqrt(error_fluid(3)^2 + error_solid(3)^2);
    error_hat(i) = sqrt(error_hat_fluid(1)^2 + error_hat_solid(1)^2);
end

r_sigma = compute_eoc(error_sigma, h_list);
r_u = compute_eoc(error_u, h_list);
r_gamma = compute_eoc(error_gamma, h_list);
r_hat = compute_eoc(error_hat, h_list);

e_o_c = [r_sigma, r_u, r_gamma, r_hat]
mean(e_o_c)
