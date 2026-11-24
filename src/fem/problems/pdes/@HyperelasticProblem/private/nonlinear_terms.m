function [m, Dm] = nonlinear_terms(obj, mesh, k, l, viscous)
% Constitutive- and symmetry-row moments of the stress remainder
% N(L) = P(I + L) − Cε(L) at the discrete gradient L_h = Σ_j l_j τ_j
% (τ_j the PEERS basis), plus, when viscous is given (viscous_rate), the
% Kelvin–Voigt stress P_v(L, L̇) of viscous.model at the rate
% l̇ = rate l + l_history (viscous_stress; ∂/∂l includes the rate):
%   m   [(ζ⁻¹N(L_h), τ_i)_K ; (N(L_h), η_i)_K]       (d2peers + d2) x 1 x N
%   Dm  ∂m/∂l                                       (d2peers + d2) x d2peers x N
% with η_i = p_i [0 1; -1 0] the rotation test functions; ζ⁻¹ is
% self-adjoint, so (ζ⁻¹N, τ) = (N, ζ⁻¹τ). l is d2peers x 1 x N. The
% volume rule is exact for degree 4k + 5, which covers the cubic Saint
% Venant–Kirchhoff remainder against the test functions. Dm is computed
% only when requested.
if nargin < 5, viscous = struct.empty; end
want_jacobian = nargout > 1;
N_ELEMENTS = size(mesh.elements, 2);
d2 = nchoosek(k + 2, 2);
d2peers = size(l, 1);
n_rows = d2peers + d2;

[lambda, mu] = LinearElasticityProblem.get_lame_parameters( ...
    obj.Data.young_modulus, obj.Data.poisson_ratio);
c1 = 1/(2*mu);
c2 = lambda/(4*mu*(lambda + mu));
material = obj.Data.material;

[xv, wv] = get_2d_reference_quadrature(2*k + 2);
r = xv(1, :).'; s = xv(2, :).';
n_pts = numel(wv);
p = dubiner_basis(k, r, s);

m = zeros(n_rows, 1, N_ELEMENTS);
if want_jacobian
    Dm = zeros(n_rows, d2peers, N_ELEMENTS);
end
chunk = 512;
for first = 1:chunk:N_ELEMENTS
    K = first:min(first + chunk - 1, N_ELEMENTS);
    nK = numel(K);
    P = PEERS_basis(k, r, s, get_affine_map(mesh, K, inverse=true));
    P = reshape(P, n_pts, d2peers, 4, nK);              % components 11, 21, 12, 22

    % test tensors [ζ⁻¹τ_i, η_i], point fastest then component
    trP = P(:, :, 1, :) + P(:, :, 4, :);
    test = zeros(n_pts, n_rows, 4, nK);
    test(:, 1:d2peers, :, :) = c1*P - c2*trP.*reshape([1 0 0 1], 1, 1, 4);
    test(:, d2peers + 1:end, 2, :) = repmat(-p, 1, 1, 1, nK);
    test(:, d2peers + 1:end, 3, :) = repmat(p, 1, 1, 1, nK);
    test = reshape(permute(test, [1 3 2 4]), n_pts*4, n_rows, nK);

    G = zeros(n_pts, 4, nK);
    for a = 1:4
        G(:, a, :) = pagemtimes(reshape(P(:, :, a, :), n_pts, d2peers, nK), l(:, :, K));
    end
    G = reshape(permute(G, [1 3 2]), n_pts*nK, 4);
    [P_G, dP_G] = HyperelasticProblem.first_piola(G, lambda, mu, material);
    [L_G, dL_G] = linear_stress(G, lambda, mu);
    if ~isempty(viscous)
        G_dot = zeros(n_pts, 4, nK);
        for a = 1:4
            G_dot(:, a, :) = pagemtimes(reshape(P(:, :, a, :), n_pts, d2peers, nK), viscous.l_history(:, :, K));
        end
        G_dot = viscous.rate*G + reshape(permute(G_dot, [1 3 2]), n_pts*nK, 4);
        [P_v, dP_v] = HyperelasticProblem.viscous_stress(G, G_dot, viscous.viscosity, ...
            viscous.rate, viscous.model);
        P_G = P_G + P_v;
        dP_G = dP_G + dP_v;
    end

    W = wv(:) .* reshape(mesh.volume(K), 1, 1, nK);   % n_pts x 1 x nK
    WN = W .* permute(reshape(P_G - L_G, n_pts, nK, 4), [1 3 2]);
    m(:, :, K) = pagemtimes(test, 'transpose', reshape(WN, n_pts*4, 1, nK), 'none');

    if want_jacobian
        D = reshape(dP_G - dL_G, n_pts, nK, 4, 4);
        DP = zeros(n_pts, 4, d2peers, nK);
        for a = 1:4
            acc = zeros(n_pts, d2peers, nK);
            for b = 1:4
                acc = acc + reshape(D(:, :, a, b), n_pts, 1, nK) .* reshape(P(:, :, b, :), n_pts, d2peers, nK);
            end
            DP(:, a, :, :) = reshape(W .* acc, n_pts, 1, d2peers, nK);
        end
        Dm(:, :, K) = pagemtimes(test, 'transpose', reshape(DP, n_pts*4, d2peers, nK), 'none');
    end
end
end

function [L, dL] = linear_stress(G, lambda, mu)
% Cε(G) = λ tr G I + μ (G + Gᵀ) and its (constant) derivative, n x 4 layout
trG = G(:, 1) + G(:, 4);
L = [lambda*trG + 2*mu*G(:, 1), mu*(G(:, 2) + G(:, 3)), ...
     mu*(G(:, 2) + G(:, 3)), lambda*trG + 2*mu*G(:, 4)];
D = [lambda + 2*mu, 0, 0, lambda;
     0, mu, mu, 0;
     0, mu, mu, 0;
     lambda, 0, 0, lambda + 2*mu];
dL = repmat(reshape(D, 1, 4, 4), size(G, 1), 1, 1);
end
