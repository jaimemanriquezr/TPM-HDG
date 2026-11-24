function [local_operators, lhs_entries, rhs_entries] = get_matrices(pde, mesh, polynomial_degree)
%-------------------------------------------------------------------------
%Description: Newton step of the gradient-based hyperelastic HDG scheme
%             at the iterate H⁰ = Data.iterate.gradient (zero without an
%             iterate), written for the new iterate (H, u, û). With
%             B(H) = (P, ∇w) − <P n, w>, C(H) = <P n, μ> (stress_terms) and
%             DB, DC their derivatives, the local problem is
%               [Mt   t_U] [H]   [t_u]       [0                          ]
%               [DB   V_U] [u] = [V_u] û  +  [V_f + DB H⁰ − B(H⁰)        ]
%             (Mt the tensor mass matrix, t_U = (u, ∇·E), t_u = <û, E n>,
%             V_U = <τ u, w> + b₀(u, w) + Newmark mass, V_u = <τ û, w>,
%             V_f the source and Newmark history), and the trace rows are
%               Σ_K DC H − <τ u, μ> + <τ û, μ> = Σ_K DC H⁰ − C(H⁰)   (+ Neumann).
%             τ is mesh.tau. Returns {A1, A2, Af} and the condensed system,
%             taken from Data.newton_operators (newton_operators) when it
%             holds them for H⁰ and mesh.tau.
%-------------------------------------------------------------------------
arguments (Input)
    pde               (1,1)   HyperelasticGradientProblem
    mesh              SimplicialGrid
    polynomial_degree (1,1)   double = pde.Order
end
k = polynomial_degree;
data = pde.Data;
N_ELEMENTS = size(mesh.elements, 2);
d1 = k + 1;
d2 = nchoosek(k + 2, 2);
H0 = zeros(4*d2, N_ELEMENTS);
if ~isempty(data.iterate), H0 = data.iterate.gradient; end
if isfield(data, "newton_operators") && ~isempty(data.newton_operators) ...
        && isequal(data.newton_operators.gradient, H0) && isequal(data.newton_operators.tau, mesh.tau)
    local_operators = data.newton_operators.local_operators;
    lhs_entries = data.newton_operators.lhs_entries;
    rhs_entries = data.newton_operators.rhs_entries;
    return
end

M = mass_matrix(mesh, k);
[Cx, Cy] = convection_matrix(mesh, k);
CxT = permute(Cx, [2 1 3]);
CyT = permute(Cy, [2 1 3]);
[tPP, nxLP, nyLP, tLP, tLL] = boundary_matrix(mesh, k);
nxPL = permute(nxLP, [2 1 3]);
nyPL = permute(nyLP, [2 1 3]);
tPL = permute(tLP, [2 1 3]);
O = zeros(d2, d2, N_ELEMENTS);
o = zeros(d2, 3*d1, N_ELEMENTS);
Oh = zeros(3*d1, d2, N_ELEMENTS);

% (H, E) + (u, ∇·E) − <û, E n> = 0
Mt = [M, O, O, O; O, M, O, O; O, O, M, O; O, O, O, M];
t_U = [CxT, O; O, CxT; CyT, O; O, CyT];
t_u = [nxPL, o; o, nxPL; nyPL, o; o, nyPL];

% momentum: B(H) + <τ(u − û), w> + b₀(u, w) (+ inertia) = (f, w)
V_U = [tPP, O; O, tPP] + data.hooke_modulus*[M, O; O, M];
V_u = [tPL, o; o, tPL];
f = data.source_term;
V_f = [test_elements(mesh, f{1}, k); test_elements(mesh, f{2}, k)];
if data.density > 0
    % Newmark-β in displacement form, as in LinearElasticityProblem
    beta = data.newmark_beta;
    dt = data.time_step;
    mass = data.density*[M, O; O, M];
    h = history(data, zeros(2*d2, N_ELEMENTS));
    V_U = V_U + mass/(beta*dt^2);
    V_f = V_f + pagemtimes(mass, reshape(h.displacement + dt*h.velocity, 2*d2, 1, []))/(beta*dt^2) ...
        + (1/2 - beta)/beta*pagemtimes(mass, reshape(h.acceleration, 2*d2, 1, []));
end

[B, C, DB, DC] = stress_terms(pde, mesh, k, H0);
H0 = reshape(H0, 4*d2, 1, N_ELEMENTS);

A1 = [Mt, t_U;
      DB, V_U];
A2 = [t_u;
      V_u];
Af = [zeros(4*d2, 1, N_ELEMENTS);
      V_f + pagemtimes(DB, H0) - B];
local_operators = {A1, A2, Af};

A3 = [DC, [-tLP, Oh; Oh, -tLP]];
A4 = [tLL, zeros(3*d1, 3*d1, N_ELEMENTS); zeros(3*d1, 3*d1, N_ELEMENTS), tLL];
c = pagemtimes(DC, H0) - C;
X = pagemldivide(A1, cat(2, A2, Af));
lhs_entries = pagemtimes(A3, X(:, 1:end-1, :)) + A4;
rhs_entries = -pagemtimes(A3, X(:, end, :)) + c;
end

function h = history(data, template)
% u⁰, v⁰, a⁰ from data.state.solid ([] or 0 = zero)
names = ["displacement", "velocity", "acceleration"];
source = [];
if isfield(data, "state") && ~isempty(data.state), source = data.state.solid; end
for i = 1:numel(names)
    value = [];
    if ~isempty(source) && isfield(source, names(i)), value = source.(names(i)); end
    if isempty(value) || isequal(value, 0), value = zeros(size(template)); end
    h.(names(i)) = value;
end
end
