function [lhs, flux, rhs] = StokesMatrices(k, mesh, nu, f, beta)
% HDG Stokes / Oseen (velocity-gradient form) as IndexedMatrix blocks over a
% heterogeneous global layout:
%   local  spaces (ids 1-4): L = Pk.to_tensor() (4 d2), u = Pk.to_vector()
%       (2 d2), p = Pk (d2), xi = DubinerBasis(0) (local mean multiplier)
%   global spaces (ids 5-7): u_hat = LegendreBasis(k).to_vector(2)
%       (skeleton), rho = DubinerBasis(0) (element pressure means),
%       mu = ConstantBasis (zero-mean multiplier)
%
% Element blocks are those of @StokesProblem/get_matrices; the hybrid
% columns/rows are reordered from component-major [x(3d1) y(3d1)] to
% per-face [x(d1); y(d1)] with the face index in dim 7 (hybrid trials) /
% dim 6 (hybrid tests), as in DarcyMatrices. The A2 (right-hand) blocks
% enter the LHS negated.
% A convective velocity beta ({bx; by} handles or Dubiner dofs) makes it
% the Oseen system of Solano & Vargas (2022, eq. 3.4):
% -(u x beta, grad v), <(uhat x beta) n, v>, -(uhat x beta) n in the
% flux, S = nu*(max|beta.n|/(2 nu) + 1).
arguments
    k double = 1
    mesh SimplicialGrid = UnitSquareMesh();
    nu = 1.0;
    f = {@(x, y) 0*(x.*y); @(x, y) 0*(x.*y)};
    beta cell = {}
end
is_oseen = ~isempty(beta);
S = nu;
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
if is_oseen
    [Cb, BbPL, BbLL, max_bn] = convective_matrix(mesh, k, beta, skip_faces = mesh.neufaces);
    S = nu*(max_bn/(2*nu) + 1);
    BbLP = permute(BbPL, [2 1 3]);
end

O  = zeros(d2, d2, N_ELEMENTS);
o  = zeros(d2, 3*d1, N_ELEMENTS);

%% local problem (tests: L rows, u rows, p(q) rows, xi(e) row)
G_L = [M, O, O, O; O, M, O, O; O, O, M, O; O, O, O, M];
G_U = [CxT, O; O, CxT; CyT, O; O, CyT];
G_u = [nxPL, o; o, nxPL; nyPL, o; o, nyPL];

V_L = -nu*[Cx, O, Cy, O; O, Cx, O, Cy];
V_U = [S*tPP, O; O, S*tPP];
V_p = [Cx; Cy];
V_u = [S*tPL, o; o, S*tPL];
if is_oseen
    V_U = V_U - [Cb, O; O, Cb];
    V_u = V_u - [BbLP, o; o, BbLP];
end
F_V = [test_elements(mesh, f{1}, k); test_elements(mesh, f{2}, k)];

q_U  = [CxT, CyT];
q_xi = M(:, 1, :);
q_u  = [nxPL, nyPL];

e_P   = M(1, :, :);
e_rho = M(1, 1, :);

lhs = IndexedArray(pg(G_L), [1 1]) + IndexedArray(pg(G_U), [2 1]) ...
    + IndexedArray(hybrid_cols(-G_u, d1), [5 1]) ...
    + IndexedArray(pg(V_L), [1 2]) + IndexedArray(pg(V_U), [2 2]) ...
    + IndexedArray(pg(V_p), [3 2]) ...
    + IndexedArray(hybrid_cols(-V_u, d1), [5 2]) ...
    + IndexedArray(pg(q_U), [2 3]) + IndexedArray(pg(q_xi), [4 3]) ...
    + IndexedArray(hybrid_cols(-q_u, d1), [5 3]) ...
    + IndexedArray(pg(e_P), [3 4]) + IndexedArray(pg(-e_rho), [6 4]);

rhs = IndexedMatrix(IndexedArray(pg(F_V), [-1 2]));

%% global problem (tests: u_hat flux rows, rho row, mu constraint row)
O3 = zeros(3*d1, d2, N_ELEMENTS);
o3 = zeros(3*d1, 1, N_ELEMENTS);
A3_L_x = [nu*nxLP, O3, nu*nyLP, O3];  A3_L_y = [O3, nu*nxLP, O3, nu*nyLP];
A3_u_x = [-S*tLP, O3];                A3_u_y = [O3, -S*tLP];
A3_p_x = -nxLP;                       A3_p_y = -nyLP;

vt = zeros(2*d1, 2*d1, 1, 1, N_ELEMENTS, 3);
for i = 1:3
    fi = (i-1)*d1 + (1:d1);
    vt(1:d1, 1:d1, 1, 1, :, i)         = reshape(S*tLL(fi, fi, :), d1, d1, 1, 1, N_ELEMENTS);
    vt(d1+(1:d1), d1+(1:d1), 1, 1, :, i) = reshape(S*tLL(fi, fi, :), d1, d1, 1, 1, N_ELEMENTS);
    if is_oseen
        vt(1:d1, 1:d1, 1, 1, :, i)         = vt(1:d1, 1:d1, 1, 1, :, i)         - reshape(BbLL(fi, fi, :), d1, d1, 1, 1, N_ELEMENTS);
        vt(d1+(1:d1), d1+(1:d1), 1, 1, :, i) = vt(d1+(1:d1), d1+(1:d1), 1, 1, :, i) - reshape(BbLL(fi, fi, :), d1, d1, 1, 1, N_ELEMENTS);
    end
end

flux = IndexedArray(hybrid_rows(A3_L_x, A3_L_y, d1), [1 5]) ...
     + IndexedArray(hybrid_rows(A3_u_x, A3_u_y, d1), [2 5]) ...
     + IndexedArray(hybrid_rows(A3_p_x, A3_p_y, d1), [3 5]) ...
     + IndexedArray(hybrid_rows(o3, o3, d1), [4 5]) ...
     + IndexedArray(vt, [5 5]) ...
     + IndexedArray(hybrid_cols([nxPL(1,:,:), nyPL(1,:,:)], d1), [5 6]) ...
     + IndexedArray(reshape(mesh.volume, 1, 1, 1, 1, N_ELEMENTS), [7 6]) ...
     + IndexedArray(reshape(mesh.volume, 1, 1, 1, 1, N_ELEMENTS), [6 7]);
end
