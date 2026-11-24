function state = get_kinematic_state(obj, solution, state)
% Advances the time-stepping state by the solved step.
%
%   state = get_kinematic_state(obj, solution, state)
%
% takes the state the step was solved from (obj.Data.state; see
% initial_state) and returns it advanced: state.time += Δt,
% state.step += 1, and state.solid updated by the time scheme.
%   newmark: with the history read from state.solid,
%       a¹ = (u¹ - u⁰ - Δt v⁰)/(β Δt²) - (1/2-β)/β a⁰
%       v¹ = v⁰ + Δt[(1-γ) a⁰ + γ a¹]
%     applied to the volume displacement (unknowns.field(2)), the
%     skeleton trace (unknowns.hybrid(1)) and, since they are linear in
%     the displacement, the stress and rotation dofs (unknowns.field(1),
%     field(3)): displacement/velocity/acceleration{,_skeleton} and
%     stress/rotation{,_velocity,_acceleration}.
%   onestep: displacement_previous ← displacement ← u¹,
%     displacement_skeleton ← û¹, stress, rotation ← the solved ones;
%     the other fields are left as they are.
arguments
    obj      (1,1) LinearElasticityProblem
    solution (1,1) FESolution
    state          struct
end
dt = obj.Data.time_step;
s  = state.solid;

u1_vol  = solution.Unknowns.field(2).dof;
u1_skel = solution.Unknowns.hybrid(1).dof;
sigma1  = solution.Unknowns.field(1).dof;
gamma1  = solution.Unknowns.field(3).dof;
if iscell(gamma1), gamma1 = cat(1, gamma1{:}); end

if obj.Data.time_scheme == "newmark"
    beta  = obj.Data.newmark_beta;
    gamma = obj.Data.newmark_gamma;

    u0_vol = default_zeros(field_of(s, "displacement"), u1_vol);
    v0_vol = default_zeros(field_of(s, "velocity"),     u1_vol);
    a0_vol = default_zeros(field_of(s, "acceleration"), u1_vol);
    a1_vol = (u1_vol - u0_vol - dt*v0_vol)/(beta*dt^2) - (1/2 - beta)/beta * a0_vol;
    v1_vol = v0_vol + dt*((1 - gamma)*a0_vol + gamma*a1_vol);

    u0_skel = default_zeros(field_of(s, "displacement_skeleton"), u1_skel);
    v0_skel = default_zeros(field_of(s, "velocity_skeleton"),     u1_skel);
    a0_skel = default_zeros(field_of(s, "acceleration_skeleton"), u1_skel);
    a1_skel = (u1_skel - u0_skel - dt*v0_skel)/(beta*dt^2) - (1/2 - beta)/beta * a0_skel;
    v1_skel = v0_skel + dt*((1 - gamma)*a0_skel + gamma*a1_skel);

    sigma_u = default_zeros(field_of(s, "stress"),                sigma1);
    sigma_v = default_zeros(field_of(s, "stress_velocity"),       sigma1);
    sigma_a = default_zeros(field_of(s, "stress_acceleration"),   sigma1);
    gamma_u = default_zeros(field_of(s, "rotation"),              gamma1);
    gamma_v = default_zeros(field_of(s, "rotation_velocity"),     gamma1);
    gamma_a = default_zeros(field_of(s, "rotation_acceleration"), gamma1);
    sigma_a1 = (sigma1 - sigma_u - dt*sigma_v)/(beta*dt^2) - (1/2 - beta)/beta * sigma_a;
    sigma_v1 = sigma_v + dt*((1 - gamma)*sigma_a + gamma*sigma_a1);
    gamma_a1 = (gamma1 - gamma_u - dt*gamma_v)/(beta*dt^2) - (1/2 - beta)/beta * gamma_a;
    gamma_v1 = gamma_v + dt*((1 - gamma)*gamma_a + gamma*gamma_a1);

    s.displacement          = u1_vol;
    s.velocity              = v1_vol;
    s.acceleration          = a1_vol;
    s.displacement_skeleton = u1_skel;
    s.velocity_skeleton     = v1_skel;
    s.acceleration_skeleton = a1_skel;
    s.stress                = sigma1;
    s.stress_velocity       = sigma_v1;
    s.stress_acceleration   = sigma_a1;
    s.rotation              = gamma1;
    s.rotation_velocity     = gamma_v1;
    s.rotation_acceleration = gamma_a1;
else
    s.displacement_previous = default_zeros(field_of(s, "displacement"), u1_vol);
    s.displacement          = u1_vol;
    s.displacement_skeleton = u1_skel;
    s.stress                = sigma1;
    s.rotation              = gamma1;
end

state.solid = s;
state.time  = state.time + dt;
state.step  = state.step + 1;
end

function value = field_of(s, name)
% a state.solid field, [] when absent
value = [];
if isfield(s, name), value = s.(name); end
end

function A = default_zeros(A, template)
    % Empty or scalar-0 fields expand to a zero array shaped like the template.
    if isempty(A) || isequal(A, 0.0)
        A = zeros(size(template));
    end
end
