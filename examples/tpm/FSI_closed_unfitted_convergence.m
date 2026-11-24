% Convergence sweep of steady FSI on a closed (circular) interface with an
% unfitted inner subdomain: both meshes are carved from one background
% grid, and the inner (solid) side keeps the elements Sigma cuts,
% protruding across the interface. The transfer paths bridge the mismatch
% (TPM). Plots the errors per degree k = 1..3 against h^(k+1), then a
% summary across degrees.
clc; clear;
for k = 1:3
    [err{k}, h{k}] = run_tests(k);
end
% Summary across degrees, in its own figure (each run_tests call opened a
% per-degree one) and against the h the errors were computed on.
figure();
loglog(h{1}, err{1}, '*-', 'DisplayName', 'P_1'); hold on
loglog(h{2}, err{2}, '*-', 'DisplayName', 'P_2')
loglog(h{3}, err{3}, '*-', 'DisplayName', 'P_3')
grid on
legend show

function [err, h_list] = run_tests(k, max_iterations, h0)
    arguments
        k
        max_iterations = 5;
        h0 = .25;
    end

    viscosity = 0.4;
    young_modulus = 1.0;
    poisson_ratio = 0.2;

    uu_fluid = "sin(pi*x)*sin(pi*y)"; %u_fluid = curl(uu_fluid)
    p_fluid = "10*(x^2 + y^2)";
    u_solid = ["sin(pi*x)*cos(pi*y)"; "cos(pi*x)*sin(pi*y)"];

    interface = CircularCurve([0.5; 0.5], 0.2, mapping="closest");

    [fluid_problem, solid_problem, interface_sources, problem_data] = ...
        FluidStructureCoupling.load_problems(k, interface, 0, 1e-3, ...
            uu_fluid, p_fluid, u_solid, ...
            0, viscosity, 0, young_modulus, poisson_ratio, 0);

    [error_sigma, error_u, error_gamma, h_list, n_elements] = deal(zeros(max_iterations, 1));
    for i = 1:max_iterations
        disp("Iteration " + i + "/" + max_iterations);
        h1 = h0 * (1/2)^(i - 1);
        disp("Generating mesh...")
        % The inner (solid) side keeps the cut elements, protruding across
        % Sigma; the crossed background admits enough face directions for
        % the staircase to follow the interface.
        [outer_mesh, inner_mesh] = StructuredSimplicialGrid(1/h1, 1/h1, "crossed") ...
            .split(interface, cut_elements="mesh_2");
        disp("Mesh generated.")

        [fluid_solution, solid_solution] = FluidStructureCoupling.solve( ...
            fluid_problem, solid_problem, outer_mesh, inner_mesh, interface, ...
            interface_sources{:});

        error_fluid = fluid_solution.compute_l2_error({problem_data.sigma_F, problem_data.u_F, problem_data.gamma_F}, {problem_data.u_F});
        error_solid = solid_solution.compute_l2_error({problem_data.sigma_S, problem_data.u_S, problem_data.gamma_S}, {problem_data.u_S});

        n_elements(i) = size(outer_mesh.elements, 2) + size(inner_mesh.elements, 2);
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

    plot_ref = @(p) loglog(h_list, error_u(1) / h_list(1)^(p) * h_list.^(p), 'k--', 'DisplayName', sprintf("$h^{%i}$", p));
    plot_ref(k+1)

    legend(ax, 'show', 'Interpreter', 'latex', 'Location', 'best', 'FontSize', 12)
    xlabel(ax, "$h$", 'Interpreter', 'latex')
    ylabel(ax, "Error", 'Interpreter', 'latex')
    title(ax, sprintf('Unfitted Stokes/Elasticity --- $\\mathcal{P}_{%i}$', k), 'Interpreter', 'latex', 'FontSize', 16);

    err = error_sigma + error_u + error_gamma;
end
