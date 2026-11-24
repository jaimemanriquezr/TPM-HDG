function [lhs, flux, rhs, geometry] = IncompressibleElasticityMatrices(k, mesh, mu, f, options)
% Incompressible elasticity (stress-velocity PEERS Stokes) as
% IndexedMatrix blocks over a heterogeneous global layout:
%   local  (ids 1-4): sigma = PEERSBasis(k), u = Pk.to_vector(),
%       gamma = Pk (scalar rotation), xi = DubinerBasis(0) (local trace
%       multiplier)
%   global (ids 5-...): u_hat = LegendreBasis(k).to_vector(2) (skeleton),
%       rho = DubinerBasis(0) (element means of tr sigma), and, under
%       pressure_mean = "constrained" only, mu_c = ConstantBasis (the
%       zero-mean multiplier; "pinned" omits it)
%
% Element blocks are those of @IncompressibleElasticityProblem/get_matrices;
% A2 blocks are negated into the LHS. With density > 0 the theta-scheme
% momentum mass terms are added (theta = 1 backward Euler, theta = 1/2
% Crank-Nicolson written as BE with dt -> theta*dt plus the
% ((1-theta)/theta) M a0 history). The histories u0 = state.velocity,
% a0 = state.acceleration (the fluid sub-struct of
% @IncompressibleElasticityProblem's state) are volume dof arrays
% (2*d2 x n_el; [] or the scalar 0 = zero).
arguments
    k double = 1
    mesh SimplicialGrid = UnitSquareMesh();
    mu = 1e-3;
    f = {@(x, y) 0*(x.*y); @(x, y) 0*(x.*y)};
    options.pressure_mean (1,1) string ...
        {mustBeMember(options.pressure_mean, ["constrained", "pinned"])} = "constrained"
    options.density   (1,1) double = 0.0
    options.time_step (1,1) double = 1.0
    options.theta     (1,1) double = 1.0
    options.state                struct = struct.empty
    % Navier-Stokes/Oseen convective velocity {bx; by} (handles or Dubiner
    % dofs), as in @IncompressibleElasticityProblem/get_matrices: kept out
    % of the stress, rho_c = density (transient) or 1 (steady)
    options.convective_velocity cell = {}
    % faces carrying the do-nothing outflow for the convective flux, in
    % addition to mesh.neufaces: faces whose componentwise boundary spec
    % has a "neumann" component must be passed here, since this builder
    % never sees the specs
    options.natural_faces (:,1) double = []
    % false leaves S = mu (the mesh tau untouched): the caller has put
    % tau(beta) on the mesh so that transmission rows built from mesh.tau
    % (InterfaceProblem) match the local problem
    options.convective_stabilization (1,1) logical = true
    % "conservative": the convective term is div(u x beta) (Solano-Vargas
    % 2022, eq. 3.4); "advective": (beta.grad) u, i.e. the same blocks
    % minus rho_c (div(beta) u, v)_K, for a beta that is not
    % divergence-free (an ALE mesh velocity); Dubiner dof arrays only
    options.convective_form (1,1) string ...
        {mustBeMember(options.convective_form, ["conservative", "advective"])} = "conservative"
    % iterate-independent element blocks of a previous call on the same
    % mesh (its fourth output), unscaled by S: M_PEERS, M_tr_PEERS,
    % M_A_PEERS, Cx, Cy, tPP, nxLP, nyLP, tLP, tLL, and M when density > 0,
    % with a key (k, n_el, n_faces, tau = mesh.tau) that must match the
    % call. Empty computes them; the fourth output returns them
    options.Cache struct = struct.empty
end
N_ELEMENTS = size(mesh.elements, 2);
d1 = k + 1;
d2 = nchoosek(k + 2, 2);
if k == 0, d2star = 0; else, d2star = k + 1; end

c1 = 1/(2*mu);
c2 = 1/(4*mu);

key = struct("k", k, "n_el", N_ELEMENTS, "n_faces", size(mesh.faces, 2), "tau", mesh.tau);
geometry = options.Cache;
if isempty(geometry)
    geometry = struct("key", key);
    [geometry.M_PEERS, geometry.M_tr_PEERS, geometry.M_A_PEERS] = PEERS_matrix(mesh, k);
    [geometry.Cx, geometry.Cy] = convection_matrix(mesh, k);
    [geometry.tPP, geometry.nxLP, geometry.nyLP, geometry.tLP, geometry.tLL] = boundary_matrix(mesh, k);
elseif ~isfield(geometry, "key") || ~isequal(geometry.key, key)
    error("IncompressibleElasticityMatrices:cacheKey", ...
        "The Cache was built for another mesh or degree (k, n_el, n_faces or tau differ).");
end
[M_PEERS, M_tr_PEERS, M_A_PEERS] = deal(geometry.M_PEERS, geometry.M_tr_PEERS, geometry.M_A_PEERS);
[Cx, Cy] = deal(geometry.Cx, geometry.Cy);
CxT = permute(Cx, [2 1 3]); CyT = permute(Cy, [2 1 3]);
[tPP, nxLP, nyLP, tLP, tLL] = deal(geometry.tPP, geometry.nxLP, geometry.nyLP, geometry.tLP, geometry.tLL);
nxPL = permute(nxLP, [2 1 3]);
nyPL = permute(nyLP, [2 1 3]);
tPL  = permute(tLP, [2 1 3]);
is_oseen = ~isempty(options.convective_velocity);
S = mu;
if is_oseen
    rho_c = 1;
    if options.density > 0; rho_c = options.density; end
    [Cb, BbPL, BbLL, max_bn] = convective_matrix(mesh, k, options.convective_velocity, skip_faces = unique([mesh.neufaces(:); options.natural_faces(:)]));
    Cb = rho_c*Cb;  BbPL = rho_c*BbPL;  BbLL = rho_c*BbLL;
    if options.convective_stabilization
        S = mu + rho_c*max_bn/2;
    end
    BbLP = permute(BbPL, [2 1 3]);
    if options.convective_form == "advective"
        Cb = Cb + rho_c*divergence_mass(mesh, k, options.convective_velocity);
    end
end

%% local problem (tests: sigma(t) rows, u(V) rows, gamma(e) row, xi(z) row)
O = zeros(d2, d2, N_ELEMENTS);
o = zeros(d2star, d2, N_ELEMENTS);
t_s = c1*M_PEERS - c2*M_tr_PEERS;
t_U = [CxT, O; O, CxT; CyT, O; O, CyT; o, o];
t_g = M_A_PEERS;
t_xi = M_tr_PEERS(:, 1, :);
Ob = zeros(d2, 3*d1, N_ELEMENTS);
ob = zeros(d2star, 3*d1, N_ELEMENTS);
t_u = [nxPL, Ob; Ob, nxPL; nyPL, Ob; Ob, nyPL; ob, ob];

V_s = -permute(t_U, [2 1 3]);
V_U = [S*tPP, O; O, S*tPP];
V_u = [S*tPL, Ob; Ob, S*tPL];
if is_oseen
    V_U = V_U - [Cb, O; O, Cb];
    V_u = V_u - [BbLP, Ob; Ob, BbLP];
end
V_f = [test_elements(mesh, f{1}, k); test_elements(mesh, f{2}, k)];

if options.density > 0
    if ~isfield(geometry, "M"), geometry.M = mass_matrix(mesh, k); end
    M = geometry.M;
    mass = options.density*[M, O; O, M];
    theta = options.theta;
    dt = options.time_step;
    h = fluid_history(options);
    V_U = V_U + mass/(theta*dt);
    if ~is_zero(h.velocity)
        u_0 = reshape(h.velocity, 2*d2, 1, []);
        V_f = V_f + pagemtimes(mass, u_0)/(theta*dt);
    end
    if theta < 1 && ~is_zero(h.acceleration)
        a_0 = reshape(h.acceleration, 2*d2, 1, []);
        V_f = V_f + (1 - theta)/theta * pagemtimes(mass, a_0);
    end
end

e_s = permute(t_g, [2 1 3]);

z_s   = permute(t_xi, [2 1 3]);
z_rho = M_PEERS(1, 1, :);

lhs = IndexedArray(pg(t_s), [1 1]) + IndexedArray(pg(t_U), [2 1]) ...
    + IndexedArray(pg(t_g), [3 1]) + IndexedArray(pg(t_xi), [4 1]) ...
    + IndexedArray(hybrid_cols(-t_u, d1), [5 1]) ...
    + IndexedArray(pg(V_s), [1 2]) + IndexedArray(pg(V_U), [2 2]) ...
    + IndexedArray(hybrid_cols(-V_u, d1), [5 2]) ...
    + IndexedArray(pg(e_s), [1 3]) ...
    + IndexedArray(pg(z_s), [1 4]) + IndexedArray(pg(-z_rho), [6 4]);

rhs = IndexedMatrix(IndexedArray(pg(V_f), [-1 2]));

%% global problem (tests: u_hat flux rows, rho row, optional mu row)
O3 = zeros(3*d1, d2, N_ELEMENTS);
o3 = zeros(3*d1, d2star, N_ELEMENTS);
A3_s_x = [nxLP, O3, nyLP, O3, o3];  A3_s_y = [O3, nxLP, O3, nyLP, o3];
A3_u_x = [-S*tLP, O3];              A3_u_y = [O3, -S*tLP];

vt = zeros(2*d1, 2*d1, 1, 1, N_ELEMENTS, 3);
for i = 1:3
    fi = (i-1)*d1 + (1:d1);
    vt(1:d1, 1:d1, 1, 1, :, i)           = reshape(S*tLL(fi, fi, :), d1, d1, 1, 1, N_ELEMENTS);
    vt(d1+(1:d1), d1+(1:d1), 1, 1, :, i) = reshape(S*tLL(fi, fi, :), d1, d1, 1, 1, N_ELEMENTS);
    if is_oseen
        vt(1:d1, 1:d1, 1, 1, :, i)           = vt(1:d1, 1:d1, 1, 1, :, i)           - reshape(BbLL(fi, fi, :), d1, d1, 1, 1, N_ELEMENTS);
        vt(d1+(1:d1), d1+(1:d1), 1, 1, :, i) = vt(d1+(1:d1), d1+(1:d1), 1, 1, :, i) - reshape(BbLL(fi, fi, :), d1, d1, 1, 1, N_ELEMENTS);
    end
end

flux = IndexedArray(hybrid_rows(A3_s_x, A3_s_y, d1), [1 5]) ...
     + IndexedArray(hybrid_rows(A3_u_x, A3_u_y, d1), [2 5]) ...
     + IndexedArray(vt, [5 5]) ...
     + IndexedArray(hybrid_cols([nxPL(1,:,:), nyPL(1,:,:)], d1), [5 6]);
if options.pressure_mean == "constrained"
    flux = flux ...
         + IndexedArray(reshape(mesh.volume, 1, 1, 1, 1, N_ELEMENTS), [7 6]) ...
         + IndexedArray(reshape(mesh.volume, 1, 1, 1, 1, N_ELEMENTS), [6 7]);
end
end

function h = fluid_history(options)
% The volume histories u0, a0 from options.state ([] when absent).
names = ["velocity", "acceleration"];
for i = 1:numel(names)
    h.(names(i)) = [];
    if isfield(options.state, names(i)) && isnumeric(options.state.(names(i)))
        h.(names(i)) = options.state.(names(i));
    end
end
end

function z = is_zero(A)
% [] and the scalar 0 stand for a zero history
z = isempty(A) || isequal(A, 0.0);
end

function D = divergence_mass(mesh, k, beta)
% (div(beta) P_b, P_a)_K, d2 x d2 x n_el, for beta = {bx; by} Dubiner dofs
assert(~isa(beta{1}, "function_handle"), "IncompressibleElasticityMatrices:advectiveHandle", ...
    "convective_form = ""advective"" needs beta as Dubiner dof arrays.");
N = size(mesh.elements, 2);
kb = find(arrayfun(@(d) nchoosek(d + 2, 2), 0:k+2) == size(beta{1}, 1)) - 1;
[rs, w] = get_2d_reference_quadrature(max(ceil((kb - 1 + 2*k)/2), 1));
p = dubiner_basis(k, rs);
[~, pr, ps] = dubiner_basis(kb, rs);
G = get_affine_map(mesh, inverse = true);
g = @(i, j) reshape(G(i, j, :), 1, N);
% d/dx = G11 d/dr + G21 d/ds, d/dy = G12 d/dr + G22 d/ds
div_beta = (pr*beta{1}).*g(1, 1) + (ps*beta{1}).*g(2, 1) ...
         + (pr*beta{2}).*g(1, 2) + (ps*beta{2}).*g(2, 2);          % nq x N
D = pagemtimes(p, 'transpose', reshape(w.*div_beta, [], 1, N).*p, 'none') ...
    .* reshape(mesh.volume, 1, 1, []);
end
