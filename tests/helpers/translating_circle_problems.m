function [fluid_problem, solid_problem, sources, problem_data, curve] = translating_circle_problems( ...
    k, curve_0, amplitude, time, time_step, uu_fluid, p_fluid, u_solid, ...
    fluid_density, fluid_viscosity, solid_density, young_modulus, poisson_ratio, rotation)
%TRANSLATING_CIRCLE_PROBLEMS  Manufactured FSI data on a translating circle.
%
%   [fp, sp, sources, pd, curve] = translating_circle_problems(k, curve_0, ...
%       A, t, dt, uu_fluid, p_fluid, u_solid, rho_F, mu, rho_S, E, nu, Theta)
%
%   Sigma(t) = Sigma_0 + d(t), d(t) = A sin(2 pi t) e_x, curve_0 a
%   CircularCurve; curve is the circle at time t. Both meshes are fixed
%   (Eulerian fluid and solid); the solid displacement u_S(x, t) is
%   measured from the translated reference, so its material derivative is
%   D_t = d_t + d'(t) d_x and the solid moves with velocity d' + D_t u_S.
%   The strings may contain t (and so d(t)).
%   With Theta (default 0) the reference also rotates about its centre,
%   phi_t(X) = c(t) + Q(theta(t)) (X - c0), c(t) = c0 + d(t) e_x,
%   theta(t) = Theta sin(2 pi t), Q the rotation matrix: Sigma(t) is the
%   same circle, and the material derivative is D_t = d_t + beta.grad with
%   the rigid velocity field beta(x, t) = d'(t) e_x + theta'(t) J (x -
%   c(t)), J = [0 -1; 1 0] (the components of u_S stay fixed-frame).
%
%   The data are load_problems' at time t on curve (fluid source, T1 with
%   the normals of Sigma(t), Dirichlet data), with the time derivatives of
%   the solid replaced by material ones, as the scheme reads them
%   (Newmark along the characteristics, LinearElasticityProblem.translate_state):
%     solid source   rho_S D_t^2 u_S - div sigma_S
%     T2             u_F - D_t u_S, the datum of the kinematic row
%                    u_F - v_S = T2, v_S the Newmark velocity D_t u_S
%                    (a physical solid gives T2 = d'(t), its rigid velocity)
%     initial_velocity, pd.u_t_S = D_t u_S;  pd.u_tt_S = D_t^2 u_S.
%   pd.T_2_t is removed: the t = 0 fluid acceleration datum of a moving
%   Sigma is not d/dt T2, so the march starts with a backward-Euler fluid
%   step (initial_state then only projects the fluid velocity).
%
%See also: FluidStructureCoupling.load_problems, LinearElasticityProblem.translate_state
d_t = sprintf("%.17g*sin(2*pi*t)", amplitude);
d_now = amplitude*sin(2*pi*time);
curve = CircularCurve(curve_0.center + [d_now; 0], curve_0.radius);

[fluid_problem, solid_problem, sources, problem_data] = FluidStructureCoupling.load_problems( ...
    k, curve, time, time_step, uu_fluid, p_fluid, u_solid, fluid_density, fluid_viscosity, ...
    solid_density, young_modulus, poisson_ratio, 0);

% material derivatives of the solid displacement
syms x y t real
u = str2sym(u_solid(:));
d_prime = diff(str2sym(d_t), t);
material = @(f) diff(f, t) + d_prime*diff(f, x);
if nargin > 13 && rotation ~= 0
    theta_prime = diff(str2sym(sprintf("%.17g*sin(2*pi*t)", rotation)), t);
    c = curve_0.center + [str2sym(d_t); 0];
    material = @(f) diff(f, t) + (d_prime - theta_prime*(y - c(2)))*diff(f, x) ...
        + theta_prime*(x - c(1))*diff(f, y);
end
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
