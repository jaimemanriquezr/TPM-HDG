function [R_u, R_hat] = gradient_operators(mesh, k, A1, A2)
% The discrete displacement gradient L_h in S_h (PEERS) of (u_h, û_h):
%   (L_h, τ)_K = −(u_h, ∇·τ)_K + <û_h, τ n>_∂K   for every τ in S_h,
% i.e. the dofs l = R_u u + R_hat û_K, from the u and û blocks of the
% constitutive row of the linear local operators A1, A2.
d2 = nchoosek(k + 2, 2);
d2peers = size(A1, 1) - 3*d2;
ix_s = 1:d2peers;
ix_u = d2peers + (1:2*d2);
M_PEERS = PEERS_matrix(mesh, k);
R_u = -pagemldivide(M_PEERS, A1(ix_s, ix_u, :));
R_hat = pagemldivide(M_PEERS, A2(ix_s, :, :));
end
