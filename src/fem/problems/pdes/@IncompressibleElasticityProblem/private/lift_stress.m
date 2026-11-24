function [sigma, gamma, div_sigma] = lift_stress(mesh, k, mu, u, u_hat)
% Viscous stress of a given velocity and trace, element by element.
%
%   [sigma, gamma, div_sigma] = lift_stress(mesh, k, mu, u, u_hat)
%
% solves the constitutive and symmetry rows of the local problem with
% (u, û) given and the trace term dropped (c1 = 1/(2μ), c2 = 0):
%   c1 (σ, τ) + (γ, τ) = <û, τ n>_∂K - (u, ∇·τ)   ∀τ
%            (σ, η) = 0                           ∀η
% so that σ = 2μ ε_h(u) (full, not deviatoric: the deviatoric operator is
% singular on traces). Linear in (u, û).
%   u          volume dofs, 2*d2 x N (x block over y block)
%   u_hat      skeleton dofs, 2*d1 x F
%   sigma      PEERS stress dofs, d2peers x N (the layout of field 1)
%   gamma      rotation dofs, d2 x N
%   div_sigma  (∇·σ, v)_K moments, 2*d2 x N
N_ELEMENTS = size(mesh.elements, 2);
d1 = k + 1;
d2 = nchoosek(k + 2, 2);
if k == 0, d2star = 0; else, d2star = k + 1; end
d2peers = 4*d2 + d2star;

[M_PEERS, ~, M_A_PEERS] = PEERS_matrix(mesh, k);
[Cx, Cy] = convection_matrix(mesh, k);
CxT = permute(Cx, [2 1 3]); CyT = permute(Cy, [2 1 3]);
[~, nxLP, nyLP] = boundary_matrix(mesh, k);
nxPL = permute(nxLP, [2 1 3]);
nyPL = permute(nyLP, [2 1 3]);

O = zeros(d2, d2, N_ELEMENTS);
o = zeros(d2star, d2, N_ELEMENTS);
t_U = [CxT, O; O, CxT; CyT, O; O, CyT; o, o];
Ob = zeros(d2, 3*d1, N_ELEMENTS);
ob = zeros(d2star, 3*d1, N_ELEMENTS);
t_u = [nxPL, Ob; Ob, nxPL; nyPL, Ob; Ob, nyPL; ob, ob];

% the trace dofs of each element's three faces, as get_local_solutions
% gathers them
ux_hat = u_hat(1:d1, :);
uy_hat = u_hat(d1+1:end, :);
u_K = reshape([reshape(ux_hat(:, mesh.boundaries), [3*d1, N_ELEMENTS]);
               reshape(uy_hat(:, mesh.boundaries), [3*d1, N_ELEMENTS])], [], 1, N_ELEMENTS);

A = [M_PEERS/(2*mu), M_A_PEERS;
     permute(M_A_PEERS, [2 1 3]), zeros(d2, d2, N_ELEMENTS)];
b = [pagemtimes(t_u, u_K) - pagemtimes(t_U, reshape(u, 2*d2, 1, N_ELEMENTS));
     zeros(d2, 1, N_ELEMENTS)];
x = pagemldivide(A, b);
sigma = reshape(x(1:d2peers, :, :), d2peers, N_ELEMENTS);
gamma = reshape(x(d2peers + 1:end, :, :), d2, N_ELEMENTS);
div_sigma = reshape(pagemtimes(t_U, "transpose", x(1:d2peers, :, :), "none"), 2*d2, N_ELEMENTS);
end
