function state = get_kinematic_state(obj, solution, state)
% Advances the Newmark state by the solved step:
%   a¹ = (u¹ − u⁰ − Δt v⁰)/(β Δt²) − (1/2 − β)/β a⁰,   v¹ = v⁰ + Δt[(1−γ) a⁰ + γ a¹]
% on the volume displacement and on the skeleton trace; gradient ← H¹.
% state.time += Δt, state.step += 1.
arguments
    obj      (1,1) HyperelasticGradientProblem
    solution (1,1) FESolution
    state          struct
end
dt = obj.Data.time_step;
beta = obj.Data.newmark_beta;
gamma = obj.Data.newmark_gamma;
s = state.solid;
u1 = solution.Unknowns.field(2).dof;
if iscell(u1), u1 = cat(1, u1{:}); end
u1_skel = solution.Unknowns.hybrid(1).dof;
if iscell(u1_skel), u1_skel = cat(1, u1_skel{:}); end
[s.displacement, s.velocity, s.acceleration] = advance(u1, s, "", dt, beta, gamma);
[s.displacement_skeleton, s.velocity_skeleton, s.acceleration_skeleton] = ...
    advance(u1_skel, s, "_skeleton", dt, beta, gamma);
s.gradient = solution.Unknowns.field(1).dof;
state.solid = s;
state.time = state.time + dt;
state.step = state.step + 1;
end

function [u1, v1, a1] = advance(u1, s, suffix, dt, beta, gamma)
get = @(name) zero_if_empty(s, name + suffix, u1);
u0 = get("displacement"); v0 = get("velocity"); a0 = get("acceleration");
a1 = (u1 - u0 - dt*v0)/(beta*dt^2) - (1/2 - beta)/beta*a0;
v1 = v0 + dt*((1 - gamma)*a0 + gamma*a1);
end

function A = zero_if_empty(s, name, template)
A = [];
if isfield(s, name), A = s.(name); end
if isempty(A) || isequal(A, 0), A = zeros(size(template)); end
end
