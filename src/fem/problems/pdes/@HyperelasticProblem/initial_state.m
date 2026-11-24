function state = initial_state(obj, mesh, k)
% Newmark state at t = obj.Data.time from the initial data.
%
%   state = initial_state(obj, mesh, k)        (k defaults to obj.Order)
%
% The state of LinearElasticityProblem.initial_state (same fields), with
% the displacement part replaced by the hyperelastic elliptic projection
% of u(0): (σ⁰, u⁰, γ⁰, û⁰) solves the static problem of this class with
% source g = −∇·P(I + ∇u(0)) + b₀ u(0) and Dirichlet datum u(0) on every
% boundary face (u(0) = Data.initial_displacement, which must evaluate on
% symbolic (x, y), as in the linear class). The velocity part (v⁰, v̂⁰ and
% its stress and rotation) is the linear elliptic projection of u̇(0), the
% consistent acceleration a⁰ is that of the momentum row (which is linear
% in this formulation: ρ (a, v) = (f, v) + (∇·σ⁰, v) − <τ(u⁰ − û⁰), v> −
% b₀ (u⁰, v), with Data.source_term = f(0) and, under "adaptive", τ
% the adaptive τ at (u⁰, û⁰)), â⁰ the face mean of its
% traces, and the stress and rotation of the acceleration are the linear
% lifts of (a⁰, â⁰). Zero initial displacement keeps the linear state.
arguments
    obj  (1,1) HyperelasticProblem
    mesh       SimplicialGrid
    k    (1,1) double = obj.Order
end
obj.Order = k;
lin = linear_problem(obj);
state = initial_state(lin, mesh, k);

w = obj.Data.initial_displacement;
syms x y real
w_sym = [sym(w{1}(x, y)); sym(w{2}(x, y))];
if all(w_sym == 0)
    return
end
[lambda, mu] = LinearElasticityProblem.get_lame_parameters( ...
    obj.Data.young_modulus, obj.Data.poisson_ratio);
P = HyperelasticProblem.first_piola_symbolic(jacobian(w_sym, [x, y]), lambda, mu, obj.Data.material);
g_sym = -[divergence(P(1, :), [x, y]); divergence(P(2, :), [x, y])] + obj.Data.hooke_modulus*w_sym;
g = cell(2, 1);
for c = 1:2
    g_c = matlabFunction(g_sym(c), Vars = [x, y]);
    g{c} = @(x, y) g_c(x, y) + 0*x;
end

static = obj;
static.Data.density = 0;
static.Data.state = struct.empty;
static.Data.iterate = struct.empty;
static.Data.source_term = g;
static.Data.boundary_conditions = boundary_condition("Dirichlet", x_data = w{1}, y_data = w{2});
mesh_D = mesh;
is_boundary = accumarray(mesh.boundaries(:), 1, [size(mesh.faces, 2), 1]) == 1;
mesh_D.faces(3, is_boundary) = mesh_D.boundary_tags("Dirichlet");
unknowns = solve(static, mesh_D, k).Unknowns;
gamma = unknowns.field(3).dof;
if iscell(gamma), gamma = cat(1, gamma{:}); end
state.solid.displacement          = unknowns.field(2).dof;
state.solid.displacement_skeleton = unknowns.hybrid(1).dof;
state.solid.stress                = unknowns.field(1).dof;
state.solid.rotation              = gamma;

% a⁰ uses the τ of the scheme at the state
mesh_tau = mesh;
if obj.Data.stabilization == "adaptive"
    mesh_tau.tau = adaptive_tau(obj, mesh, k, struct( ...
        "displacement", state.solid.displacement, ...
        "displacement_skeleton", state.solid.displacement_skeleton));
end
a0  = lin.consistent_acceleration(mesh_tau, k, state);
a0h = skeleton_mean(mesh, k, a0);
[sigma_a, gamma_a] = lin.lift_stress(mesh, k, a0, a0h);
state.solid.acceleration          = a0;
state.solid.acceleration_skeleton = a0h;
state.solid.stress_acceleration   = sigma_a;
state.solid.rotation_acceleration = gamma_a;
end
