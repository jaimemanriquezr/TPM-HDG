% Steady fluid-structure interaction on a flat interface with an explicit
% gap: the two fitted meshes retreat d from y0 and the transfer paths
% bridge the band (TPM). Prints the L2 errors on each side and plots both
% solutions.
clear;
close all;

k = 1;
h = 1/4;

y0 = .5;
d = 1e-3;

viscosity = 0.4;
young_modulus = 1.0;
poisson_ratio = 0.25;

uu_fluid = "x*y";
p_fluid = "x^3 + y^3";
u_solid = ["x"; "y"];

interface = FlatCurveHorizontal(y0);
[upper_mesh, lower_mesh] = fitted_split(interface, h, d);

[fluid_problem, solid_problem, interface_sources, problem_data] = ...
    FluidStructureCoupling.load_problems(k, interface, 0, 1e-3, ...
        uu_fluid, p_fluid, u_solid, ...
        0, viscosity, 0, young_modulus, poisson_ratio, 0);

[fluid_solution, solid_solution] = FluidStructureCoupling.solve( ...
    fluid_problem, solid_problem, upper_mesh, lower_mesh, interface, ...
    interface_sources{:}, quiet=true);

error_l2_fluid = fluid_solution.compute_l2_error({problem_data.sigma_F, problem_data.u_F, problem_data.gamma_F}, {problem_data.u_F})
error_l2_solid = solid_solution.compute_l2_error({problem_data.sigma_S, problem_data.u_S, problem_data.gamma_S}, {problem_data.u_S})

plot(fluid_solution);
plot(solid_solution);
