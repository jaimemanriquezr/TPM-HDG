function handle = biot_path_gradient(k, young_modulus, poisson_ratio, biot_alpha, kappa)
%BIOT_PATH_GRADIENT  Transfer-path gradient reconstruction for HDG Biot.
%
%   handle = biot_path_gradient(k, young_modulus, poisson_ratio, biot_alpha, kappa)
%
%   returns the function-handle form (b) of
%   TransferEssentialBoundaryConditions.path_gradient for the Biot element
%   of BiotMatrices, whose trace is (eta_hat_x, eta_hat_y, p_hat):
%
%     Gm = handle(problem, mesh, rs, K_list, m, n_local),
%     Gm (N x n_local x n_b x 3).
%
%   The local columns follow the Biot layout [sigma (dS); eta (2 d2);
%   gamma (d2); q (2 d2); p (d2)]. With the total stress sigma, the
%   compliance A tau = c1 tau - c2 tr(tau) I and c_b = alpha/(2(lambda+mu)),
%
%     grad(eta) = A(sigma + alpha p I) + gamma = A sigma + c_b p I + gamma,
%     grad(p)   = -kappa^{-1} q,
%
%   so the components of [G_phi m] are
%
%     a = 1, 2:  [(A sigma_phi + c_b p_phi I + gamma_phi) m]_a
%                (sigma, gamma and p columns; zero on eta and q),
%     a = 3:     -(kappa^{-1} q_phi) . m   (q columns),
%
%   with gamma = g [0 1; -1 0] the rotation of LinearElasticityMatrices.
%   The transferred condition is u_hat_a(x) + int_0^1 [G m]_a ds = g_a(x_bar).
%
%See also: TransferEssentialBoundaryConditions, BiotMatrices
arguments
    k             (1,1) double
    young_modulus (1,1) double
    poisson_ratio (1,1) double
    biot_alpha    (1,1) double
    kappa         (2,2) double = eye(2)
end
[lambda, mu] = LinearElasticityProblem.get_lame_parameters(young_modulus, poisson_ratio);
c1 = 1/(2*mu);
c2 = lambda/(4*mu*(lambda + mu));
c_b = biot_alpha/(2*(lambda + mu));
handle = @(problem, mesh, rs, K_list, m, n_local) ...
    biot_gradient(mesh, rs, K_list, m, n_local, k, c1, c2, c_b, kappa);
end

function Gm = biot_gradient(mesh, rs, K_list, m, n_local, k, c1, c2, c_b, kappa)
dS = PEERSBasis(k).Dimension;
d2 = nchoosek(k + 2, 2);
assert(n_local == dS + 6*d2, 'biot_path_gradient:layout', ...
    'expected %d local dofs (Biot layout, k = %d), got %d', dS + 6*d2, k, n_local);
N = size(rs, 2);  n_b = size(rs, 3);
i_gamma = dS + 2*d2 + (1:d2);
i_qx = dS + 3*d2 + (1:d2);
i_qy = dS + 4*d2 + (1:d2);
i_p = dS + 5*d2 + (1:d2);

G_maps = get_affine_map(mesh, K_list(:), inverse = true, output = "array");
a = reshape(kappa.' \ reshape(m, 2, []), 2, N, n_b);     % kappa^{-T} m

Gm = zeros(N, n_local, n_b, 3);
for b = 1:n_b
    m_b = reshape(m(:, :, b), 2, N).';                   % (N x 2)
    % compliance applied to the PEERS stress, contracted with m
    P = PEERSBasis.evaluate_basis(k, rs(1, :, b), rs(2, :, b), G_maps(:, :, b));
    trP = P(:, :, 1, 1) + P(:, :, 2, 2);
    AP = c1*P - c2*trP.*reshape(eye(2), 1, 1, 2, 2);
    for c = 1:2
        Gm(:, 1:dS, b, c) = AP(:, :, c, 1).*m_b(:, 1) + AP(:, :, c, 2).*m_b(:, 2);
    end
    phi = DubinerBasis.evaluate_basis(k, rs(1, :, b).', rs(2, :, b).');   % (N x d2)
    % rotation: gamma m = g [m2; -m1]
    Gm(:, i_gamma, b, 1) = phi.*m_b(:, 2);
    Gm(:, i_gamma, b, 2) = -phi.*m_b(:, 1);
    % pore pressure in the displacement gradient: c_b p m
    Gm(:, i_p, b, 1) = c_b*phi.*m_b(:, 1);
    Gm(:, i_p, b, 2) = c_b*phi.*m_b(:, 2);
    % Darcy: grad(p) m = -(kappa^{-T} m) . q
    Gm(:, i_qx, b, 3) = -phi.*reshape(a(1, :, b), N, 1);
    Gm(:, i_qy, b, 3) = -phi.*reshape(a(2, :, b), N, 1);
end
end
