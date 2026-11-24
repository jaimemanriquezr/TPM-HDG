function [fluid_problem, solid_problem, problem_data, reference, sources] = cavity_problems(k, time, time_step, options)
%CAVITY_PROBLEMS  Lid-driven cavity with an elastic disc, physical data.
%
%   [fp, sp, pd, reference, sources] = cavity_problems(k, t, dt, ...)
%
%   The data of lagrangian_carved_march (its problems option, as
%   @(t) cavity_problems(k, t, dt, ...)) for the disc in a lid-driven
%   cavity of Zhao et al. (2008) and Sugiyama et al. (2011): the unit
%   cavity, a Stokes fluid (density rho_f, viscosity mu_f) driven by the
%   lid velocity on y = 1,
%     u = (r(t), 0),   r(t) = (1 - cos(pi t/ramp_time))/2
%   for t < ramp_time and 1 after (a C^1 start from rest; ramp_time = 0
%   is the impulsive start r = 1), u = 0 on the other walls: the uniform
%   lid of the references, discontinuous at the top corners (the datum is
%   projected face-wise on the Dirichlet faces at Gauss points, which on
%   the side walls all lie below y = 1, so the lid faces carry r(t) and
%   the side faces 0). lid_profile = "regularized" replaces 1 by
%   16 x^2 (1 - x)^2, which vanishes at the corners. The disc of radius R about center is the reference Sigma_0,
%   a neo-Hookean HyperelasticProblem with shear modulus G
%   (E = 2 G (1 + nu)), density rho_s and the Kelvin-Voigt viscosity
%   solid_viscosity; viscous_model other than "material" is stored as the
%   solid's Data.viscous_model (left unset for "material").
%   Sources, interface data (T1, T2) and initial data are zero. pd holds
%   zero fields in the layout of lagrangian_march_problems, so that the
%   march's "errors" are the norms of the discrete solution and its
%   geometry error is the largest interface displacement.
%
%See also: lagrangian_carved_march, lagrangian_march_problems
arguments
    k (1,1) double
    time (1,1) double
    time_step (1,1) double
    options.center (2,1) double = [0.6; 0.5]
    options.radius (1,1) double = 0.2
    options.rho_f (1,1) double = 1
    options.mu_f (1,1) double = 0.01
    options.rho_s (1,1) double = 1
    options.shear_modulus (1,1) double = 0.1
    options.nu (1,1) double = 0.49
    options.solid_viscosity (1,1) double = 0.01
    options.ramp_time (1,1) double {mustBeNonnegative} = 0.1
    options.lid_profile (1,1) string {mustBeMember(options.lid_profile, ["regularized", "uniform"])} = "uniform"
    options.viscous_model (1,1) string {mustBeMember(options.viscous_model, ["material", "spatial"])} = "material"
end
reference = CircularCurve(options.center, options.radius);
zero = @(x, y) 0*x;
T_r = options.ramp_time;
if time < T_r
    ramp = (1 - cos(pi*time/T_r))/2;
    ramp_rate = pi/(2*T_r)*sin(pi*time/T_r);
else
    ramp = 1;
    ramp_rate = 0;
end
on_lid = @(y) abs(y - 1) < 1e-10;
lid = @(x, y, r) r*16*x.^2.*(1 - x).^2.*on_lid(y);
if options.lid_profile == "uniform"
    lid = @(x, y, r) r*on_lid(y) + 0*x;
end
fluid_problem = IncompressibleElasticityProblem(k, time_step = time_step, ...
    viscosity = options.mu_f, density = options.rho_f, ...
    domain_measure = 1 - pi*options.radius^2, pressure_mean = "pinned", ...
    datum_dirichlet = {@(x, y) lid(x, y, ramp); zero}, ...
    datum_dirichlet_rate = {@(x, y) lid(x, y, ramp_rate); zero});
E = 2*options.shear_modulus*(1 + options.nu);
solid_problem = HyperelasticProblem(k, material = "neo_hookean", young_modulus = E, ...
    poisson_ratio = options.nu, density = options.rho_s, time_step = time_step, time = time, ...
    solid_viscosity = options.solid_viscosity);
if options.viscous_model ~= "material"
    solid_problem.Data.viscous_model = options.viscous_model;
end
problem_data = struct("sigma_F", {{zero, zero; zero, zero}}, "u_F", {{zero; zero}}, ...
    "gamma_F", zero, "sigma_S", {{zero, zero; zero, zero}}, "u_S", {{zero; zero}}, ...
    "gamma_S", zero);
sources = @(curve) {{zero; zero}, {zero; zero}};
end
