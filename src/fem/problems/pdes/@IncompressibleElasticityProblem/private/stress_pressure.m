function p = stress_pressure(mesh, k, sigma)
% Dubiner dofs (d2 x N) of p = -tr(σ)/2, the element L² projection of the
% trace of the PEERS stress dofs sigma (d2peers x N): the Dubiner blocks
% contribute -(σ_11 + σ_22)/2, the bubbles through their trace moments.
N_ELEMENTS = size(mesh.elements, 2);
d2 = nchoosek(k + 2, 2);
[~, M_tr_PEERS] = PEERS_matrix(mesh, k);
M = mass_matrix(mesh, k);
moments = pagemtimes(M_tr_PEERS(1:d2, :, :), reshape(sigma, [], 1, N_ELEMENTS));
p = -reshape(pagemldivide(M, moments), d2, N_ELEMENTS)/2;
end
