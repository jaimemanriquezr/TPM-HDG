function [B, C, DB, DC] = stress_terms(obj, mesh, k, H)
% Element contributions of the Piola stress P(I + H_h):
%   B(w) = (P, ∇w)_K − <P n, w>_∂K      momentum rows,  2*d2 x 1 x N
%   C(μ) = <P n, μ>_∂K                  trace rows,     6*d1 x 1 x N
%   DB, DC = ∂B/∂H, ∂C/∂H               2*d2 (6*d1) x 4*d2 x N (when requested)
% H holds the dofs of H_h (4*d2 x N, components 11, 21, 12, 22). The
% rules are exact for degree 4k + 5 (volume and faces).
want_jacobian = nargout > 2;
N_ELEMENTS = size(mesh.elements, 2);
d1 = k + 1;
d2 = nchoosek(k + 2, 2);
nH = 4*d2;
[lambda, mu] = LinearElasticityProblem.get_lame_parameters( ...
    obj.Data.young_modulus, obj.Data.poisson_ratio);
material = obj.Data.material;

% reference points: the volume rule, then the three faces (local face e
% is (0, 1-s), (1-s, 0), (1-s, s), as in boundary_matrix)
[xv, wv] = get_2d_reference_quadrature(2*k + 2);
[s, ws] = get_1d_reference_quadrature(2*k + 2);
s = s(:); ws = ws(:);
nv = numel(wv); nf = numel(s);
r_all = [xv(1, :).'; 0*s; 1 - s; 1 - s];
s_all = [xv(2, :).'; 1 - s; 0*s; s];
n_pts = nv + 3*nf;
face_rows = @(e) nv + (e-1)*nf + (1:nf);

[p, pr, ps] = dubiner_basis(k, r_all, s_all);          % n_pts x d2
pr(nv + 1:end, :) = 0; ps(nv + 1:end, :) = 0;          % gradients at volume points only
pf = p; pf(1:nv, :) = 0;                               % values at face points only
l_face = {legendre_basis(k, s), legendre_basis(k, 1 - s)};
oriented_normals = mesh.normals .* kron(mesh.orientation, [1; 1]);

B = zeros(2*d2, 1, N_ELEMENTS);
C = zeros(6*d1, 1, N_ELEMENTS);
if want_jacobian
    DB = zeros(2*d2, nH, N_ELEMENTS);
    DC = zeros(6*d1, nH, N_ELEMENTS);
end

chunk = 1024;
for first = 1:chunk:N_ELEMENTS
    K = first:min(first + chunk - 1, N_ELEMENTS);
    nK = numel(K);
    Ginv = get_affine_map(mesh, K, inverse=true);
    W = reshape([wv(:) .* reshape(mesh.volume(K), 1, nK); repmat(ws, 3, nK)], n_pts, 1, nK);
    nx = zeros(n_pts, 1, nK); ny = zeros(n_pts, 1, nK);
    for e = 1:3
        nx(face_rows(e), 1, :) = repmat(reshape(oriented_normals(2*e - 1, K), 1, 1, nK), nf, 1);
        ny(face_rows(e), 1, :) = repmat(reshape(oriented_normals(2*e, K), 1, 1, nK), nf, 1);
    end

    % momentum test tensors: row (c, i) is e_c ⊗ ∇w_i in the volume and
    % −e_c ⊗ n w_i on the faces (slots 11, 21, 12, 22)
    gwx = pr .* reshape(Ginv(1, 1, :), 1, 1, nK) + ps .* reshape(Ginv(2, 1, :), 1, 1, nK) - nx.*pf;
    gwy = pr .* reshape(Ginv(1, 2, :), 1, 1, nK) + ps .* reshape(Ginv(2, 2, :), 1, 1, nK) - ny.*pf;
    Psi_B = zeros(n_pts, 4, 2*d2, nK);
    Psi_B(:, 1, 1:d2, :) = reshape(gwx, n_pts, 1, d2, nK);
    Psi_B(:, 3, 1:d2, :) = reshape(gwy, n_pts, 1, d2, nK);
    Psi_B(:, 2, d2 + 1:end, :) = reshape(gwx, n_pts, 1, d2, nK);
    Psi_B(:, 4, d2 + 1:end, :) = reshape(gwy, n_pts, 1, d2, nK);
    Psi_B = reshape(Psi_B, n_pts*4, 2*d2, nK);

    % trace test tensors: row (c, e, i) is e_c ⊗ n μ_i on face e
    perm = mesh.permutation(:, K);
    lf = zeros(n_pts, 3*d1, nK);
    for e = 1:3
        for o = 1:2
            is_o = find(perm(e, :) == o);
            lf(face_rows(e), (e-1)*d1 + (1:d1), is_o) = repmat(l_face{o}, 1, 1, numel(is_o));
        end
    end
    Psi_C = zeros(n_pts, 4, 6*d1, nK);
    Psi_C(:, 1, 1:3*d1, :) = reshape(nx.*lf, n_pts, 1, 3*d1, nK);
    Psi_C(:, 3, 1:3*d1, :) = reshape(ny.*lf, n_pts, 1, 3*d1, nK);
    Psi_C(:, 2, 3*d1 + 1:end, :) = reshape(nx.*lf, n_pts, 1, 3*d1, nK);
    Psi_C(:, 4, 3*d1 + 1:end, :) = reshape(ny.*lf, n_pts, 1, 3*d1, nK);
    Psi_C = reshape(Psi_C, n_pts*4, 6*d1, nK);

    % H_h and P at the points
    HK = reshape(H(:, K), d2, 4, nK);
    G = pagemtimes(p, HK);                               % n_pts x 4 x nK
    G = reshape(permute(G, [1 3 2]), n_pts*nK, 4);
    [P, dP] = HyperelasticGradientProblem.first_piola(G, lambda, mu, material);
    WP = W .* permute(reshape(P, n_pts, nK, 4), [1 3 2]);
    WP = reshape(WP, n_pts*4, 1, nK);
    B(:, :, K) = pagemtimes(Psi_B, 'transpose', WP, 'none');
    C(:, :, K) = pagemtimes(Psi_C, 'transpose', WP, 'none');

    if want_jacobian
        % ∂P_a/∂H_(b,j) = D_ab p_j
        D = reshape(dP, n_pts, nK, 4, 4);
        DP = zeros(n_pts, 4, d2, 4, nK);
        for a = 1:4
            for b = 1:4
                DP(:, a, :, b, :) = reshape(W .* reshape(D(:, :, a, b), n_pts, 1, nK) .* p, ...
                    n_pts, 1, d2, 1, nK);
            end
        end
        DP = reshape(DP, n_pts*4, nH, nK);
        DB(:, :, K) = pagemtimes(Psi_B, 'transpose', DP, 'none');
        DC(:, :, K) = pagemtimes(Psi_C, 'transpose', DP, 'none');
    end
end
end
