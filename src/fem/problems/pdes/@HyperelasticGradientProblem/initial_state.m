function state = initial_state(obj, mesh, k)
% Newmark state at t = obj.Data.time from the initial data.
%
%   state = initial_state(obj, mesh, k)        (k defaults to obj.Order)
%
% (H⁰, u⁰, û⁰) is the hyperelastic elliptic projection of u(0): the
% static solve of this class with source −∇·P(I + ∇u(0)) + b₀ u(0) and
% Dirichlet datum u(0) on every boundary face; (v⁰, v̂⁰) the same
% projection of u̇(0) with the linear law. The consistent acceleration
% solves ρ (a⁰, w) = (f, w) − B(H⁰) − <τ(u⁰ − û⁰), w> − b₀ (u⁰, w)
% (Data.source_term = f(0)), â⁰ is the face mean of its traces. The
% handles must evaluate on symbolic (x, y). Needs density > 0.
% Returns state.time, state.step = 0, state.solid (fields as in the
% class help), state.fluid = state.diagnostics = struct.empty.
arguments
    obj  (1,1) HyperelasticGradientProblem
    mesh       SimplicialGrid
    k    (1,1) double = obj.Order
end
obj.Order = k;
data = obj.Data;
N_ELEMENTS = size(mesh.elements, 2);
N_FACES = size(mesh.faces, 2);
d1 = k + 1;
d2 = nchoosek(k + 2, 2);
[H0, u0, u0h] = projection(obj, mesh, k, data.initial_displacement, data.material);
[~, v0, v0h] = projection(obj, mesh, k, data.initial_velocity, "linear");
solid = struct("displacement", u0, "velocity", v0, "acceleration", zeros(2*d2, N_ELEMENTS), ...
    "displacement_skeleton", u0h, "velocity_skeleton", v0h, ...
    "acceleration_skeleton", zeros(2*d1, N_FACES), "gradient", H0);
state = struct("time", data.time, "step", 0, "solid", solid, ...
    "fluid", struct.empty, "diagnostics", struct.empty);
assert(data.density > 0, "HyperelasticGradientProblem:density", ...
    "initial_state requires density > 0.");
% consistent acceleration: ρ M a⁰ = −(momentum residual of the static problem)
static = obj;
static.Data.density = 0;
static.Data.state = struct.empty;
static.Data.iterate = struct("gradient", H0);
static_mesh = mesh;
static_mesh.tau = stabilization_tau(static, mesh, k, H0);
local_operators = get_matrices(static, static_mesh, k);
[A1, A2, Af] = deal(local_operators{:});
ix_u = 4*d2 + (1:2*d2);
X = reshape([H0; u0], [], 1, N_ELEMENTS);
r = pagemtimes(A1(ix_u, :, :), X) - pagemtimes(A2(ix_u, :, :), element_traces(mesh, u0h, d1)) - Af(ix_u, :, :);
M = mass_matrix(mesh, k);
O = zeros(d2, d2, N_ELEMENTS);
a0 = reshape(pagemldivide(data.density*[M, O; O, M], -r), 2*d2, N_ELEMENTS);
state.solid.acceleration = a0;
state.solid.acceleration_skeleton = skeleton_mean(mesh, k, a0);
end

function [H, u, u_hat] = projection(obj, mesh, k, w, material)
% the static solve of material with the source of w and Dirichlet data w
syms x y real
w_sym = [sym(w{1}(x, y)); sym(w{2}(x, y))];
d2 = nchoosek(k + 2, 2);
if all(w_sym == 0)
    H = zeros(4*d2, size(mesh.elements, 2));
    u = zeros(2*d2, size(mesh.elements, 2));
    u_hat = zeros(2*(k + 1), size(mesh.faces, 2));
    return
end
[lambda, mu] = LinearElasticityProblem.get_lame_parameters(obj.Data.young_modulus, obj.Data.poisson_ratio);
G = jacobian(w_sym, [x, y]);
if material == "linear"
    P = lambda*trace(G)*eye(2) + mu*(G + G.');
else
    P = HyperelasticProblem.first_piola_symbolic(G, lambda, mu, material);
end
g_sym = -[divergence(P(1, :), [x, y]); divergence(P(2, :), [x, y])] + obj.Data.hooke_modulus*w_sym;
g = cell(2, 1);
for c = 1:2
    g_c = matlabFunction(g_sym(c), Vars = [x, y]);
    g{c} = @(x, y) g_c(x, y) + 0*x;
end
static = obj;
static.Data.material = material;
static.Data.density = 0;
static.Data.state = struct.empty;
static.Data.iterate = struct.empty;
static.Data.source_term = g;
static.Data.boundary_conditions = boundary_condition("Dirichlet", x_data = w{1}, y_data = w{2});
mesh_D = mesh;
is_boundary = accumarray(mesh.boundaries(:), 1, [size(mesh.faces, 2), 1]) == 1;
mesh_D.faces(3, is_boundary) = mesh_D.boundary_tags("Dirichlet");
unknowns = solve(static, mesh_D, k).Unknowns;
H = unknowns.field(1).dof;
u = unknowns.field(2).dof;
if iscell(u), u = cat(1, u{:}); end
u_hat = unknowns.hybrid(1).dof;
if iscell(u_hat), u_hat = cat(1, u_hat{:}); end
end
