function [lhs, flux, rhs] = AccelerationMatrices(k, mesh, rho, F, options)
% Velocity-hybridized HDG for the acceleration problem at a fixed time
% (the pressure Poisson / projection problem of the stress-velocity fluid),
%   rho a + grad p = F,   div a = 0   in Omega,
% as IndexedMatrix blocks over a heterogeneous global layout:
%   local  (ids 1-3): a = Pk.to_vector() (2 d2), p = Pk (d2),
%       xi = DubinerBasis(0) (local mean multiplier)
%   global (ids 4-...): a_hat = LegendreBasis(k).to_vector(2) (skeleton),
%       rho_K = DubinerBasis(0) (element pressure means), and, under
%       pressure_mean = "constrained" only, mu_c = ConstantBasis (the
%       zero-mean multiplier; "pinned" omits it)
%
% The blocks are those of StokesMatrices with the gradient unknown and its
% rows deleted and the momentum mass rho (a, v)_K added:
%   rho (a, v)_K + (grad p, v)_K + <S (a - a_hat), v>_dK = (F, v)_K
%   (a, grad q)_K - <a_hat.n, q>_dK + (xi, q)_K          = 0
%   (p, 1)_K = rho_K |K|
% and the numerical flux rows, in total-flux form with a given stress
% sigma_v (the viscous part of the stress at that time),
%   sum_K <-p n - S (a - a_hat), mu>_dK = -sum_K <sigma_v n, mu>_dK,
% so that sigma_v n - p n - S (a - a_hat) is the normal trace of the
% total stress; on traction faces the boundary datum adds <q_N, mu>. A
% further known normal flux g (e.g. a convective flux) adds
% sum_K <g, mu>_dK to the right-hand side of the flux rows.
%
%   F           (F, v)_K moments, 2*d2 x n_el (x block over y block); [] = 0
%   S           stabilization*(mesh tau), stabilization defaulting to rho
%               (S then has the units of a traction per velocity)
%   viscous_stress  PEERS dofs of sigma_v, d2peers x n_el (the layout of
%               the stress-velocity fluid's field 1; its bubbles have zero
%               normal trace); [] = 0
%   flux_moments    <g, mu>_e moments per element and local face e,
%               6*d1 x n_el, component-major ([x(3 d1); y(3 d1)], the
%               row layout of the low-level A3 blocks); [] = 0
arguments
    k double
    mesh SimplicialGrid
    rho (1,1) double
    F double = []
    options.pressure_mean (1,1) string ...
        {mustBeMember(options.pressure_mean, ["constrained", "pinned"])} = "constrained"
    options.stabilization (1,1) double = rho
    options.viscous_stress double = []
    options.flux_moments double = []
end
S = options.stabilization;
N_ELEMENTS = size(mesh.elements, 2);
d1 = k + 1;
d2 = nchoosek(k + 2, 2);

M = mass_matrix(mesh, k);
[Cx, Cy] = convection_matrix(mesh, k);
CxT = permute(Cx, [2 1 3]); CyT = permute(Cy, [2 1 3]);
[tPP, nxLP, nyLP, tLP, tLL] = boundary_matrix(mesh, k);
nxPL = permute(nxLP, [2 1 3]);
nyPL = permute(nyLP, [2 1 3]);
tPL  = permute(tLP, [2 1 3]);

O  = zeros(d2, d2, N_ELEMENTS);
o  = zeros(d2, 3*d1, N_ELEMENTS);

%% local problem (tests: a(v) rows, p(q) rows, xi(e) row)
V_U = rho*[M, O; O, M] + [S*tPP, O; O, S*tPP];
V_p = [Cx; Cy];
V_u = [S*tPL, o; o, S*tPL];

q_U  = [CxT, CyT];
q_xi = M(:, 1, :);
q_u  = [nxPL, nyPL];

e_P   = M(1, :, :);
e_rho = M(1, 1, :);

lhs = IndexedArray(pg(V_U), [1 1]) + IndexedArray(pg(V_p), [2 1]) ...
    + IndexedArray(hybrid_cols(-V_u, d1), [4 1]) ...
    + IndexedArray(pg(q_U), [1 2]) + IndexedArray(pg(q_xi), [3 2]) ...
    + IndexedArray(hybrid_cols(-q_u, d1), [4 2]) ...
    + IndexedArray(pg(e_P), [2 3]) + IndexedArray(pg(-e_rho), [5 3]);

if isempty(F)
    F = zeros(2*d2, N_ELEMENTS);
end
rhs = IndexedMatrix(IndexedArray(reshape(F, 2*d2, 1, 1, 1, N_ELEMENTS), [-1 1]));

%% global problem (tests: a_hat flux rows, rho row, optional mu row)
O3 = zeros(3*d1, d2, N_ELEMENTS);
o3 = zeros(3*d1, 1, N_ELEMENTS);
A3_u_x = [-S*tLP, O3];                A3_u_y = [O3, -S*tLP];
A3_p_x = -nxLP;                       A3_p_y = -nyLP;

vt = zeros(2*d1, 2*d1, 1, 1, N_ELEMENTS, 3);
for i = 1:3
    fi = (i-1)*d1 + (1:d1);
    vt(1:d1, 1:d1, 1, 1, :, i)           = reshape(S*tLL(fi, fi, :), d1, d1, 1, 1, N_ELEMENTS);
    vt(d1+(1:d1), d1+(1:d1), 1, 1, :, i) = reshape(S*tLL(fi, fi, :), d1, d1, 1, 1, N_ELEMENTS);
end

flux = IndexedArray(hybrid_rows(A3_u_x, A3_u_y, d1), [1 4]) ...
     + IndexedArray(hybrid_rows(A3_p_x, A3_p_y, d1), [2 4]) ...
     + IndexedArray(hybrid_rows(o3, o3, d1), [3 4]) ...
     + IndexedArray(vt, [4 4]) ...
     + IndexedArray(hybrid_cols([nxPL(1,:,:), nyPL(1,:,:)], d1), [4 5]);
if options.pressure_mean == "constrained"
    flux = flux ...
         + IndexedArray(reshape(mesh.volume, 1, 1, 1, 1, N_ELEMENTS), [6 5]) ...
         + IndexedArray(reshape(mesh.volume, 1, 1, 1, 1, N_ELEMENTS), [5 6]);
end

% -sum_K <sigma_v n, mu>_dK: (sigma_v n)_x = sigma_11 n_x + sigma_12 n_y,
% (sigma_v n)_y = sigma_21 n_x + sigma_22 n_y (PEERS blocks 1, 3 and 2, 4)
if ~isempty(options.viscous_stress)
    sigma = reshape(options.viscous_stress, [], 1, N_ELEMENTS);
    blk = @(b) sigma((b-1)*d2 + (1:d2), :, :);
    t_x = pagemtimes(nxLP, blk(1)) + pagemtimes(nyLP, blk(3));
    t_y = pagemtimes(nxLP, blk(2)) + pagemtimes(nyLP, blk(4));
    rhs = rhs + IndexedArray(hybrid_rows(-t_x, -t_y, d1), [-1 4]);
end
if ~isempty(options.flux_moments)
    g = reshape(options.flux_moments, 6*d1, 1, N_ELEMENTS);
    rhs = rhs + IndexedArray(hybrid_rows(g(1:3*d1, :, :), g(3*d1+1:end, :, :), d1), [-1 4]);
end
end
