function tau = adaptive_tau(obj, mesh, k, iterate)
% Face-wise adaptive stabilization (Terrana, Nguyen, Bonet, Peraire,
% CMAME 352 (2019), eq. 18):
%   τ_F = ρ_max(∂P/∂F) / (2 L_c),
% ρ_max the largest eigenvalue of the tangent ∂P/∂F (symmetric for a
% hyperelastic law) over the quadrature points of face F, seen from every
% element sharing F, at F = I + L_h of the iterate (u, û); L_c =
% Data.characteristic_length. Without an iterate L_h = 0 and τ is the
% linear value 2(μ + λ)/(2 L_c) (plane strain). Returns the 3 x N array
% of mesh.tau (equal on both sides of a face).
arguments
    obj  (1,1) HyperelasticProblem
    mesh       SimplicialGrid
    k    (1,1) double
    iterate    struct
end
[lambda, mu] = LinearElasticityProblem.get_lame_parameters( ...
    obj.Data.young_modulus, obj.Data.poisson_ratio);
L_c = obj.Data.characteristic_length;
N_ELEMENTS = size(mesh.elements, 2);
N_FACES = size(mesh.faces, 2);
if isempty(iterate)
    tau = (2*mu + 2*lambda)/(2*L_c) * ones(3, N_ELEMENTS);
    return
end
d1 = k + 1;
d2 = nchoosek(k + 2, 2);

% the discrete gradient of the iterate
operators = linear_operators(obj, mesh, k);
R_u = operators.R_u;
R_hat = operators.R_hat;
d2peers = size(R_u, 1);
u = reshape(iterate.displacement, 2*d2, 1, N_ELEMENTS);
u_K = element_traces(mesh, iterate.displacement_skeleton, d1);
l = pagemtimes(R_u, u) + pagemtimes(R_hat, u_K);

% the three faces of the reference element, (0, 1-s), (1-s, 0), (1-s, s)
s = get_1d_reference_quadrature(2*k + 2);
s = s(:); nf = numel(s);
r_all = [0*s; 1 - s; 1 - s];
s_all = [1 - s; 0*s; s];
n_pts = 3*nf;

rho = zeros(3, N_ELEMENTS);
chunk = 512;
for first = 1:chunk:N_ELEMENTS
    K = first:min(first + chunk - 1, N_ELEMENTS);
    nK = numel(K);
    P = reshape(PEERS_basis(k, r_all, s_all, get_affine_map(mesh, K, inverse=true)), ...
        n_pts, d2peers, 4, nK);
    G = zeros(n_pts, 4, nK);
    for a = 1:4
        G(:, a, :) = pagemtimes(reshape(P(:, :, a, :), n_pts, d2peers, nK), l(:, :, K));
    end
    G = reshape(permute(G, [1 3 2]), n_pts*nK, 4);
    [~, dP] = HyperelasticProblem.first_piola(G, lambda, mu, obj.Data.material);
    D = permute(dP, [2 3 1]);                         % 4 x 4 x points
    D = (D + pagetranspose(D))/2;
    rho_points = reshape(max(pageeig(D), [], 1), nf, 3, nK);
    rho(:, K) = reshape(max(rho_points, [], 1), 3, nK);
end
% the largest value over the elements sharing each face
rho_face = accumarray(mesh.boundaries(:), rho(:), [N_FACES, 1], @max);
tau = reshape(rho_face(mesh.boundaries), 3, N_ELEMENTS)/(2*L_c);
end
