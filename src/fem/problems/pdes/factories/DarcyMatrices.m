function [lhs, flux] = DarcyMatrices(k, mesh, kappa)
% HDG Darcy as IndexedMatrix blocks over a heterogeneous global layout:
%   local  spaces (ids 1-2): u = Pk.to_vector() (2 d2), flux;
%       p = Pk (d2), pressure
%   global space  (id 3):    p_hat = LegendreBasis(k) (skeleton), pressure
%       trace
%
% A block's id is [trial test]. lhs holds the blocks tested by the local
% spaces (the element equations, test ids 1-2); flux holds the blocks
% tested by the trace (test id 3), the numerical-flux row. kappa is the
% constant 2-by-2 diffusion tensor (default eye(2)); the stabilization is
% mesh.tau. The hybrid blocks carry the face index in dim 7 (trials) /
% dim 6 (tests).
    arguments
        k double
        mesh SimplicialGrid
        kappa = eye(2);
    end
    t = permute(mesh.tau, [3:6 2 1]);
    n = permute(mesh.exterior_unit_normals, [5 6 1 2 4 3]);

    P_k = DubinerBasis(k);
    L_k = LegendreBasis(k);

    mass = pages(mass_matrix(mesh, P_k, P_k));
    M = pagekron(inv(kappa), mass);

    cv = convection_matrix(mesh, P_k, P_k);
    conv = cat(3, pages(cv.x), pages(cv.y));
    C = reshape(conv, size(conv, 1), [], 1, 1, size(conv, 5));

    PP = face_pages(boundary_matrix(mesh, P_k, P_k));
    LP = face_pages(boundary_matrix(mesh, P_k, L_k));
    n_LP = n .* LP;
    N_LP = reshape(n_LP, size(n_LP, 1), [], 1, 1, size(n_LP, 5), 3);

    ll = pages(mass_matrix(mesh, L_k, L_k));
    LL = reshape(ll(:,:,:,:,mesh.boundaries.'), size(ll,1),size(ll,2),1,1,[],3);

    k_q_r_dx = IndexedArray(M, [1 1]);
    u_div_r_dx = IndexedArray(pagetranspose(C), [2 1]);
    N_PL = permute(pagetranspose(N_LP), [1:5 7 6]);
    U_r_n_dS = IndexedArray(N_PL, [3 1]);

    div_q_w_dx = IndexedArray(C, [1 2]);
    tau_u_w_dS = IndexedArray(sum(t .* PP, 6), [2 2]);
    tPL = permute(t .* pagetranspose(LP), [1:5 7 6]);
    tau_U_w_dS = IndexedArray(tPL, [3 2]);

    q_n_V_dS = IndexedArray(N_LP, [1 3]);
    tau_u_V_dS = IndexedArray(t .* LP, [2 3]);
    tau_U_V_dS = IndexedArray(t .* LL, [3 3]);

    lhs = k_q_r_dx - u_div_r_dx + U_r_n_dS ...
        + div_q_w_dx + tau_u_w_dS - tau_U_w_dS;
    flux = -(q_n_V_dS + tau_u_V_dS - tau_U_V_dS);
end
