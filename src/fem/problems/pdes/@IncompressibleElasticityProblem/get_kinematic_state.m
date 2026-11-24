function state = get_kinematic_state(obj, solution, state, options)
% Advances the time-stepping state by the solved theta-scheme step.
%
%   state = get_kinematic_state(obj, solution, state)
%
% takes the state the step was solved from (obj.Data.state; see
% initial_state; an empty state is the zero history) and returns it
% advanced: state.time += Δt, state.step += 1 and, with u⁰, a⁰ read from
% state.fluid,
%   velocity          = u¹                      (unknowns.field(2))
%   acceleration      = (u¹ - u⁰)/(θ Δt) - (1-θ)/θ a⁰
%   velocity_skeleton = û¹                      (unknowns.hybrid(1))
%   stress            = σ¹                      (unknowns.field(1), PEERS dofs)
%   pressure          = -tr(σ¹)/2               (Dubiner dofs, d2 x N)
% (θ = 1 reduces to the backward-Euler identity a¹ = (u¹ - u⁰)/Δt). For a
% steady problem (density = 0) only the solved fields are updated.
%
% theta defaults to the θ the matrices were assembled with (get_theta).
arguments
    obj      (1,1) IncompressibleElasticityProblem
    solution (1,1) FESolution
    state          struct
    options.theta  (1,1) double = NaN
end
transient = obj.Data.density > 0;
if transient
    theta = options.theta;
    if isnan(theta), theta = obj.get_theta(); end
    dt = obj.Data.time_step;
end
u1 = solution.Unknowns.field(2).dof;

if isempty(state)
    state = struct("time", obj.Data.time, "step", 0, "solid", struct.empty, ...
        "fluid", struct.empty, "diagnostics", struct.empty);
end
s = state.fluid;
if isempty(s), s = struct(); end
sigma1 = solution.Unknowns.field(1).dof;
if transient
    u0 = default_zeros(field_of(s, "velocity"),     u1);
    a0 = default_zeros(field_of(s, "acceleration"), u1);
    s.acceleration = (u1 - u0)/(theta*dt) - (1 - theta)/theta * a0;
end
s.velocity          = u1;
s.velocity_skeleton = solution.Unknowns.hybrid(1).dof;
s.stress            = sigma1;
s.pressure          = stress_pressure(solution.Mesh, solution.Order, sigma1);
state.fluid = s;
if transient
    state.time = state.time + dt;
    state.step = state.step + 1;
end
end

function value = field_of(s, name)
% a state.fluid field, [] when absent
value = [];
if isfield(s, name), value = s.(name); end
end

function A = default_zeros(A, template)
    % Empty or scalar-0 fields expand to a zero array shaped like the template.
    if isempty(A) || isequal(A, 0.0)
        A = zeros(size(template));
    end
end
