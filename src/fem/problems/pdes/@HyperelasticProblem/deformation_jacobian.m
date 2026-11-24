function J = deformation_jacobian(pde, mesh, k, iterate)
% det F = det(I + L_h) at the volume quadrature points of nonlinear_terms.
%
%   J = deformation_jacobian(pde, mesh, k, iterate)
%
% L_h is the discrete gradient of iterate (struct with displacement and
% displacement_skeleton, the layout of Data.iterate), the field at which
% the stress remainder and its Newton step are evaluated; min(J, [], "all")
% > 0 is the condition for first_piola of the neo-Hookean law to be
% defined there. J is n_points x N (one column per element).
arguments
    pde     (1,1) HyperelasticProblem
    mesh          SimplicialGrid
    k       (1,1) double
    iterate       struct
end
operators = linear_operators(pde, mesh, k);
N_ELEMENTS = size(mesh.elements, 2);
d1 = k + 1;
d2 = nchoosek(k + 2, 2);
d2peers = size(operators.R_u, 1);
u = reshape(iterate.displacement, 2*d2, 1, N_ELEMENTS);
u_K = element_traces(mesh, iterate.displacement_skeleton, d1);
l = pagemtimes(operators.R_u, u) + pagemtimes(operators.R_hat, u_K);

xv = get_2d_reference_quadrature(2*k + 2);
r = xv(1, :).'; s = xv(2, :).';
n_pts = numel(r);
J = zeros(n_pts, N_ELEMENTS);
chunk = 512;
for first = 1:chunk:N_ELEMENTS
    K = first:min(first + chunk - 1, N_ELEMENTS);
    nK = numel(K);
    P = reshape(PEERS_basis(k, r, s, get_affine_map(mesh, K, inverse=true)), n_pts, d2peers, 4, nK);
    G = zeros(n_pts, 4, nK);
    for a = 1:4
        G(:, a, :) = pagemtimes(reshape(P(:, :, a, :), n_pts, d2peers, nK), l(:, :, K));
    end
    J(:, K) = reshape((1 + G(:, 1, :)).*(1 + G(:, 4, :)) - G(:, 2, :).*G(:, 3, :), n_pts, nK);
end
end
