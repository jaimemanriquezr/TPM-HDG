function state = initial_state(obj, mesh, k)
% Time-stepping state at t = obj.Data.time from the initial data.
%
%   state = initial_state(obj, mesh, k)        (k defaults to obj.Order)
%
% Discretizes Data.initial_displacement, Data.initial_velocity (cells of
% handles {η_x; η_y}, {v_x; v_y}) and Data.initial_pressure (handle) by
% projection:
%   η⁰ = Π η(0), v⁰ = Π v(0), p⁰ = Π p(0) onto the elements (L2proj2D),
%   η̂⁰ = P_M η(0), p̂⁰ = P_M p(0) onto the faces (L2projSkeleton1d),
% lifts the effective stress and the rotation from (η⁰, η̂⁰) with the
% elasticity local problem of the same E, ν (LinearElasticityProblem
% lift_stress: A σ_e + γ = ∇η in the discrete sense), and stores the
% total stress σ⁰ = σ_e - α p⁰ I (BiotProblem.effective_stress with -α).
% No equation is solved: with storage s0 = 0 the initial data should
% satisfy the mass equation's constraint themselves.
%
% Returns state.time = Data.time, state.step = 0,
% state.fluid = struct.empty, state.diagnostics = struct.empty, and
% state.biot with the dof arrays (one column per element or face)
%   displacement_skeleton  2*d1 x n_faces   η̂⁰
%   stress                 dS x n_el        σ⁰ (PEERS layout)
%   rotation               d2 x n_el        γ⁰
%   pressure               d2 x n_el        p⁰
%   pressure_skeleton      d1 x n_faces     p̂⁰
%   displacement           2*d2 x n_el      η⁰ (volume; read by the inertia)
%   velocity               2*d2 x n_el      v⁰ = Π Data.initial_velocity
% and displacement_skeleton_previous, stress_previous, rotation_previous,
% pressure_previous, displacement_previous, velocity_previous = [] (no
% level n-1: a "bdf2" step from this state is an Euler step, see
% biot_step_history).
%
%See also: BiotProblem/get_kinematic_state, biot_step_history,
%          LinearElasticityProblem/lift_stress
arguments
    obj  (1,1) BiotProblem
    mesh       SimplicialGrid
    k    (1,1) double = obj.Order
end
data = obj.Data;
eta0     = L2proj2D(mesh, data.initial_displacement, k);
eta0_hat = L2projSkeleton1d(mesh, data.initial_displacement, k);
v0       = L2proj2D(mesh, data.initial_velocity, k);
p0       = L2proj2D(mesh, data.initial_pressure, k);
p0_hat   = L2projSkeleton1d(mesh, data.initial_pressure, k);

elasticity = LinearElasticityProblem(k, young_modulus = data.young_modulus, ...
    poisson_ratio = data.poisson_ratio);
[sigma_e, gamma] = elasticity.lift_stress(mesh, k, eta0, eta0_hat);
sigma = BiotProblem.effective_stress(sigma_e, p0, -data.biot_alpha, k);

biot = struct( ...
    "displacement_skeleton",          eta0_hat, ...
    "stress",                         sigma, ...
    "rotation",                       gamma, ...
    "pressure",                       p0, ...
    "pressure_skeleton",              p0_hat, ...
    "displacement",                   eta0, ...
    "velocity",                       v0, ...
    "displacement_skeleton_previous", [], ...
    "stress_previous",                [], ...
    "rotation_previous",              [], ...
    "pressure_previous",              [], ...
    "displacement_previous",          [], ...
    "velocity_previous",              []);
state = struct("time", data.time, "step", 0, "biot", biot, ...
    "fluid", struct.empty, "diagnostics", struct.empty);
end
