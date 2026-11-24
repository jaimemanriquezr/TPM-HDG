function state = initial_state(obj, mesh, k)
% Time-stepping state at t = obj.Data.time from the initial data.
%
%   state = initial_state(obj, mesh, k)        (k defaults to obj.Order)
%
% Discretizes Data.initial_displacement and Data.initial_velocity (cells
% of handles {u_x; u_y}) into volume dofs u⁰, v⁰, skeleton dofs û⁰, v̂⁰
% and the stress/rotation dofs (σ⁰, γ⁰) of the displacement.
%   onestep: the projections u⁰ = Π u(0), v⁰ = Π u̇(0) onto the elements
%     (L2proj2D) and û⁰ = P_M u(0), v̂⁰ = P_M u̇(0) onto the faces
%     (L2projSkeleton1d), with (σ⁰, γ⁰) lifted from (u⁰, û⁰)
%     (lift_stress); displacement_previous is [], which makes the first
%     step the Taylor start from u⁰, v⁰.
%   newmark: the HDG elliptic projections of u(0) and u̇(0): for w = u(0)
%     (resp. u̇(0)), (σ_h, w_h, γ_h, ŵ_h) solves the static problem of
%     this class (same E, ν, hooke_modulus b₀, stabilization and mesh)
%     with source g = -∇·Cε(w) + b₀ w and Dirichlet datum w on every
%     boundary face; (u⁰, û⁰, σ⁰, γ⁰) and (v⁰, v̂⁰) with their
%     stress_velocity, rotation_velocity are those solutions. g is derived
%     symbolically, so each handle must evaluate on symbolic (x, y) to a
%     symbolic scalar (elementwise operations only, constants written as
%     c + 0*x); otherwise the error
%     LinearElasticityProblem:initialDataNotSymbolic is raised. Then the
%     consistent acceleration a⁰ (consistent_acceleration, which reads
%     Data.source_term: set it to f(0) before the call; needs density > 0)
%     satisfies ρ a⁰ = Π (f(0) - g): when the data satisfy the equation
%     at t = 0, a⁰ = Π ü(0), where the projections of the onestep branch
%     give an O(h^k) face residual. Its skeleton â⁰ is the face L²
%     projection of the mean of the one-sided traces of a⁰ (one-sided on
%     boundary faces), and (σ, γ) of the acceleration are lifted from
%     (a⁰, â⁰).
% Returns state.time = Data.time, state.step = 0, state.solid (fields as
% in the LinearElasticityProblem help; [] = not used by the scheme),
% state.fluid = struct.empty, state.diagnostics = struct.empty.
arguments
    obj  (1,1) LinearElasticityProblem
    mesh       SimplicialGrid
    k    (1,1) double = obj.Order
end
data = obj.Data;
is_newmark = data.time_scheme == "newmark";
if is_newmark
    assert(data.density > 0, "LinearElasticityProblem:density", ...
        "initial_state with time_scheme = ""newmark"" requires density > 0.");
    [u0, u0h, sigma_u, gamma_u] = elliptic_projection(obj, mesh, k, ...
        data.initial_displacement, "initial_displacement");
    [v0, v0h, sigma_v, gamma_v] = elliptic_projection(obj, mesh, k, ...
        data.initial_velocity, "initial_velocity");
else
    u0  = L2proj2D(mesh, data.initial_displacement, k);
    u0h = L2projSkeleton1d(mesh, data.initial_displacement, k);
    v0  = L2proj2D(mesh, data.initial_velocity, k);
    v0h = L2projSkeleton1d(mesh, data.initial_velocity, k);
    [sigma_u, gamma_u] = obj.lift_stress(mesh, k, u0, u0h);
end

solid = struct( ...
    "displacement",          u0,  ...
    "velocity",              v0,  ...
    "acceleration",          [],  ...
    "displacement_previous", [],  ...
    "displacement_skeleton", u0h, ...
    "velocity_skeleton",     v0h, ...
    "acceleration_skeleton", [],  ...
    "stress",                sigma_u, ...
    "rotation",              gamma_u, ...
    "stress_velocity",       [],  ...
    "rotation_velocity",     [],  ...
    "stress_acceleration",   [],  ...
    "rotation_acceleration", []);
state = struct("time", data.time, "step", 0, "solid", solid, ...
    "fluid", struct.empty, "diagnostics", struct.empty);

if is_newmark
    a0  = obj.consistent_acceleration(mesh, k, state);
    a0h = skeleton_mean(mesh, k, a0);
    [sigma_a, gamma_a] = obj.lift_stress(mesh, k, a0, a0h);
    state.solid.acceleration          = a0;
    state.solid.acceleration_skeleton = a0h;
    state.solid.stress_velocity       = sigma_v;
    state.solid.rotation_velocity     = gamma_v;
    state.solid.stress_acceleration   = sigma_a;
    state.solid.rotation_acceleration = gamma_a;
end
end

function [w_h, w_hat, sigma, gamma] = elliptic_projection(obj, mesh, k, w, name)
% HDG elliptic projection of the vector field w = {w_x; w_y} (handles):
% the static solution of this problem with source -∇·Cε(w) + b₀ w (C the
% plane-strain Hooke law of get_matrices) and Dirichlet datum w on every
% boundary face. Returns the volume (2*d2 x N), skeleton (2*d1 x F),
% stress (PEERS layout) and rotation dofs.
syms x y real
w_sym = sym(zeros(2, 1));
for c = 1:2
    is_symbolic = false;
    try
        value = w{c}(x, y);
        is_symbolic = isa(value, "sym") && isscalar(value);
    catch
    end
    if ~is_symbolic
        error("LinearElasticityProblem:initialDataNotSymbolic", ...
            "Data.%s{%d} must evaluate on symbolic (x, y) to a symbolic scalar " + ...
            "(the newmark initial state derives its elliptic projection " + ...
            "symbolically): build it from elementwise operations, constants " + ...
            "as c + 0*x.", name, c);
    end
    w_sym(c) = value;
end
[lambda, mu] = LinearElasticityProblem.get_lame_parameters( ...
    obj.Data.young_modulus, obj.Data.poisson_ratio);
grad_w = jacobian(w_sym, [x, y]);
epsilon = (grad_w + grad_w.')/2;
sigma_sym = 2*mu*epsilon + lambda*trace(epsilon)*eye(2);
g_sym = -[divergence(sigma_sym(1, :), [x, y]); divergence(sigma_sym(2, :), [x, y])] ...
    + obj.Data.hooke_modulus*w_sym;
g = cell(2, 1);
for c = 1:2
    g_c = matlabFunction(g_sym(c), Vars = [x, y]);
    g{c} = @(x, y) g_c(x, y) + 0*x;
end

static = obj;
static.Data.density = 0;
static.Data.state = struct.empty;
static.Data.source_term = g;
static.Data.boundary_conditions = boundary_condition("Dirichlet", ...
    x_data = w{1}, y_data = w{2});
% every boundary face carries the Dirichlet tag
mesh_D = mesh;
is_boundary = accumarray(mesh.boundaries(:), 1, [size(mesh.faces, 2), 1]) == 1;
mesh_D.faces(3, is_boundary) = mesh_D.boundary_tags("Dirichlet");
unknowns = solve(static, mesh_D, k).Unknowns;
w_h   = unknowns.field(2).dof;
w_hat = unknowns.hybrid(1).dof;
sigma = unknowns.field(1).dof;
gamma = unknowns.field(3).dof;
if iscell(gamma), gamma = cat(1, gamma{:}); end
end

function a_hat = skeleton_mean(mesh, k, a)
% face L² projection of the mean of the one-sided traces of the volume
% dofs a (2*d2 x N): on each face F, Σ_K <a_K, L>_F / n_F against the face
% mass, n_F = the number of elements sharing F
d1 = k + 1;
d2 = size(a, 1)/2;
N_ELEMENTS = size(mesh.elements, 2);
N_FACES = size(mesh.faces, 2);
% <P_j, L_i>_e and <L_j, L_i>_e per element and local face e (unit weight)
[~, ~, ~, LP, LL] = boundary_matrix(add_stabilization(mesh, 1), k);
incidence = sparse(1:3*N_ELEMENTS, mesh.boundaries(:), 1, 3*N_ELEMENTS, N_FACES);
face_mass = zeros(d1, d1, 3, N_ELEMENTS);
for e = 1:3
    block_e = (e-1)*d1 + (1:d1);
    face_mass(:, :, e, :) = reshape(LL(block_e, block_e, :), d1, d1, 1, N_ELEMENTS);
end
face_mass = reshape(reshape(face_mass, d1*d1, 3*N_ELEMENTS)*incidence, d1, d1, N_FACES);
a_hat = zeros(2*d1, N_FACES);
for c = 1:2
    a_c = reshape(a((c-1)*d2 + (1:d2), :), d2, 1, N_ELEMENTS);
    moments = reshape(pagemtimes(LP, a_c), d1, 3*N_ELEMENTS)*incidence;
    a_hat((c-1)*d1 + (1:d1), :) = reshape( ...
        pagemldivide(face_mass, reshape(moments, d1, 1, N_FACES)), d1, N_FACES);
end
end
