function L_h = gradient_from_stress(mesh, k, sigma_dof, gamma_dof, mu)
arguments
    mesh SimplicialGrid
    k (1,1) double
    sigma_dof (:,:) double   % PEERS stress dofs (4*d2 + d1) x Nelts
    gamma_dof (:,:) double   % rotation dofs, d2 x Nelts
    mu (1,1) double          % shear viscosity
end
%-------------------------------------------------------------------------
%Description: Velocity gradient of the stress-velocity (PEERS) fluid,
%               grad u = (1/2mu) sigma - (1/4mu) tr(sigma) I + gamma
%             (the incompressible limit of the constitutive row, as in
%             trace_extension_coefficients), L2-projected element by
%             element onto Dubiner P_k so it can feed postprocess_velocity.
%Output:    - L_h: {L11, L12; L21, L22}, d2 x Nelts each
%-------------------------------------------------------------------------
N_ELEMENTS = size(mesh.elements, 2);
c1 = 1/(2*mu);  c2 = 1/(4*mu);
[~, rs, w] = quadpoints_elements_2d(mesh, k + 1, output="cell");
G = get_affine_map(mesh, inverse = true);
Pk = dubiner_basis(k, rs{:});                       % nq x d2
M = (w .* Pk).' * Pk;                               % reference mass (volume cancels)
A = asym_tensor_dubiner_basis(k, rs{1}.', rs{2}.'); % nq x d2 x 2 x 2

L_h = cell(2, 2);
[L_h{:}] = deal(zeros(size(Pk, 2), N_ELEMENTS));
for K = 1:N_ELEMENTS
    P = PEERS_basis(k, rs{1}.', rs{2}.', G(:, :, K));         % nq x nb x 2 x 2
    sigma = reshape(tensorprod(P, sigma_dof(:, K), 2, 1), [], 2, 2);
    gamma = reshape(tensorprod(A, gamma_dof(:, K), 2, 1), [], 2, 2);
    trace_sigma = sigma(:, 1, 1) + sigma(:, 2, 2);
    grad_u = c1*sigma + gamma;
    grad_u(:, 1, 1) = grad_u(:, 1, 1) - c2*trace_sigma;
    grad_u(:, 2, 2) = grad_u(:, 2, 2) - c2*trace_sigma;
    for i = 1:2
        for j = 1:2
            L_h{i, j}(:, K) = M \ ((w .* Pk).' * grad_u(:, i, j));
        end
    end
end
end
