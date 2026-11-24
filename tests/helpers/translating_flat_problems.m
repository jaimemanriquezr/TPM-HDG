function [fluid_problem, solid_problem, sources, problem_data, curve] = translating_flat_problems( ...
    k, y0, amplitude, time, time_step, uu_fluid, p_fluid, u_solid, ...
    fluid_density, fluid_viscosity, solid_density, young_modulus, poisson_ratio)
%TRANSLATING_FLAT_PROBLEMS  Manufactured FSI data on a vertically moving flat interface.
%
%   [fp, sp, sources, pd, curve] = translating_flat_problems(k, y0, A, t, dt, ...
%       uu_fluid, p_fluid, u_solid, rho_F, mu, rho_S, E, nu)
%
%   Sigma(t) = {y = y0 + d(t)}, d(t) = A sin(2 pi t); fluid above, solid
%   below. The flat counterpart of translating_circle_problems: the solid
%   displacement u_S(x, t) is measured from the reference translated by
%   d(t) e_y, so its material derivative is D_t = d_t + d'(t) d_y; the
%   solid source, the kinematic datum T2 = u_F - D_t u_S, the initial
%   velocity and pd.u_t_S, pd.u_tt_S use D_t (see that file).
%
%See also: translating_circle_problems, FluidStructureCoupling.load_problems
d_t = sprintf("%.17g*sin(2*pi*t)", amplitude);
curve = FlatCurveHorizontal(y0 + amplitude*sin(2*pi*time));

[fluid_problem, solid_problem, sources, problem_data] = FluidStructureCoupling.load_problems( ...
    k, curve, time, time_step, uu_fluid, p_fluid, u_solid, fluid_density, fluid_viscosity, ...
    solid_density, young_modulus, poisson_ratio, 0);

syms x y t real
u = str2sym(u_solid(:));
d_prime = diff(str2sym(d_t), t);
material = @(f) diff(f, t) + d_prime*diff(f, y);
D_u  = material(u);
D2_u = material(D_u);
at_time = @(e) arrayfun(@(c) vectorized(matlabFunction(subs(c, t, time), Vars = [x, y])), ...
    e(:), UniformOutput = false);
extra_source = at_time(solid_density*(D2_u - diff(u, t, 2)));
extra_T2 = at_time(-(D_u - diff(u, t)));
D_u  = at_time(D_u);
D2_u = at_time(D2_u);
for c = 1:2
    f = solid_problem.Data.source_term{c};  g = extra_source{c};
    solid_problem.Data.source_term{c} = @(x, y) f(x, y) + g(x, y);
    T2 = sources{2}{c};  g = extra_T2{c};
    sources{2}{c} = @(x, y) T2(x, y) + g(x, y);
end
solid_problem.Data.initial_velocity = D_u;
problem_data.f_S = solid_problem.Data.source_term;
problem_data.T_2 = sources{2};
problem_data.u_t_S = D_u;
problem_data.u_tt_S = D2_u;
problem_data = rmfield(problem_data, "T_2_t");
end

function h = vectorized(f)
% broadcasts spatially constant components to the array shape
h = @(x, y) f(x, y) + 0*x;
end
