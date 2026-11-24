% Steady fluid-structure interaction on a closed interface: the solid
% occupies the disk, the fluid the rest of the unit square, and both
% subdomains are meshed fitted to Sigma (fitted_split, closed topology).
% Prints the L2 errors on each side and plots both solutions.
clear;
close all;

k = 2;
h = 1/6;

I = CircularCurve([0.5; 0.5], 0.2, mapping="closest");

viscosity = 0.4;
young_modulus = 1.0;
poisson_ratio = 0.25;

uu_fluid = "sin(pi*x)*cos(pi*y)";
p_fluid = "0";
u_solid = ["pi*sin(pi*x)*sin(pi*y) + 1"; "pi*cos(pi*x)*cos(pi*y) - 1"];

[outer_mesh, inner_mesh] = fitted_split(I, h, 0);

[fluid_problem, solid_problem, interface_sources, problem_data] = ...
    FluidStructureCoupling.load_problems(k, I, 0, 1e-3, ...
        uu_fluid, p_fluid, u_solid, ...
        0, viscosity, 0, young_modulus, poisson_ratio, 0);

[fluid_solution, solid_solution] = FluidStructureCoupling.solve( ...
    fluid_problem, solid_problem, outer_mesh, inner_mesh, I, ...
    interface_sources{:}, quiet=true);

error_l2_fluid = fluid_solution.compute_l2_error({problem_data.sigma_F, problem_data.u_F, problem_data.gamma_F}, {problem_data.u_F})
error_l2_solid = solid_solution.compute_l2_error({problem_data.sigma_S, problem_data.u_S, problem_data.gamma_S}, {problem_data.u_S})

plot(fluid_solution);
plot(solid_solution);
