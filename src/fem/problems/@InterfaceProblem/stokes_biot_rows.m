function [Hc, Fc] = stokes_biot_rows(obj, k, S1, S2, mesh_1, mesh_2, ...
    part_1, part_2, off, n_dof, extra_base, sq, tq)
% Stokes-Biot transmission rows over the induced patch partition (see
% InterfaceProblem.stokes_biot). Side 1 = PEERS stress-velocity fluid,
% local [sigma_f(dS); u_f(2 d2); gamma_f(d2); xi(1)], global [u_hat_f,
% rho]; side 2 = Biot, local [sigma_p(dS); eta(2 d2); gamma_p(d2);
% q(2 d2); p(d2)], one global trace [eta_hat_x; eta_hat_y; p_hat] (3 d1
% per face). Local solves are x = S1f - S12 x_hat, so a functional R on
% the local unknowns enters as -R*S12 in the matrix and -R*S1f in the rhs.
%
% Per patch e, with nu_p the unit normal of Sigma out of side 2, nu_f =
% -nu_p, t = (-nu_p2, nu_p1), tau_i = part_i.tau(e), tau_s =
% stabilization_solid * tau_2, tau_p = stabilization_darcy * ||kappa|| *
% tau_2 (as in BiotMatrices), m_i the
% transfer path vector of side i:
%   Eh(sigma_f) = E(sigma_f) nu_f - tau_1 mu_f (u_f - u_hat_f)|face
%   Eh(sigma_p) = E(sigma_p) nu_p - tau_s (eta - eta_hat)|face
%   Eh(q_hat)   = E(q).nu_p + tau_p (p - p_hat)|face
%   Th(u_hat_f) = u_hat_f + int_path (A_f E(sigma_f) + E(gamma_f)) m_1
%   Th(eta_hat) = eta_hat + int_path (A E(sigma_p) + c_b E(p) I + E(gamma_p)) m_2
%   Tp(p_hat)   = p_hat - int_path (kappa^-1 E(q)) . m_2
% and the rows
%   R_f on the fluid trace dofs (vector Legendre test, side-1 measure):
%       <Eh(sigma_f) + Eh(sigma_p), mu> = <T_1, mu>
%   R_n on the eta_hat_x dofs (scalar Legendre test, side-2 measure):
%       <nu_p.(cdt Th(u_hat_f) - Th(eta_hat)) - cdt Eh(q_hat), l>
%           = <-nu_p.h* + cdt nu_p.T_2, l>
%   R_t on the eta_hat_y dofs:
%       <t.Eh(sigma_f) + a t.(cdt Th(u_hat_f) - Th(eta_hat)), l>
%           = <-a t.h* + t.T_4, l>,   a = friction_bjs / cdt
%   R_p on the p_hat dofs:
%       <Tp(p_hat) + nu_f.Eh(sigma_f), l> = <T_3, l>
% h* = eta_hat* + Ext(eta*) is the transported displacement history:
% the coupling's history_skeleton gathered per patch plus its transport
% onto Sigma (the coupling's history_extension, or
% trace_extension_coefficients of the effective stress sigma* + alpha p* I
% and rotation gamma* on the Biot side of this partition; none when those
% histories are zero). Extrapolated terms are evaluated at the patch
% points on Sigma (quadrature.interface_2d), face terms on the face
% (quadrature.local_2d / local).
% The pore-pressure level scalar beta sits at extra_base + 1: column
% -cdt <beta, l> in R_n, and the row
%   |Omega_P,h|^-1 (p, 1)_{Omega_P,h} = pore_pressure_mean
% on the recovered Biot pressure (see InterfaceProblem.stokes_biot;
% InterfaceProblem.solve eliminates beta when the level is fixed).
d1 = k + 1;
d2 = nchoosek(k + 2, 2);
dS = PEERSBasis(k).Dimension;
c = obj.coupling;
mu_f = c.viscosity;
c1_1 = 1/(2*mu_f);
c2_1 = 1/(4*mu_f);
[lambda_p, mu_p] = LinearElasticityProblem.get_lame_parameters( ...
    c.young_modulus, c.poisson_ratio);
c1_2 = 1/(2*mu_p);
c2_2 = lambda_p/(2*mu_p*(2*lambda_p + 2*mu_p));
c_b = c.biot_alpha/(2*(lambda_p + mu_p));
cdt = c.effective_time_step;
a_dof = extra_base + 1;
a_bjs = c.friction_bjs/cdt;
kappa = c.diffusion;
tau_p_factor = BiotProblem.darcy_tau_factor(c.stabilization_darcy, kappa);
[~, wS]  = lgwt(sq, 0, 1);
[~, wTP] = lgwt(tq, 0, 1);
n_p = size(part_1.faces, 3);

tr1 = obj.problems(1).GlobalSpaces([obj.problems(1).GlobalSpaces.Support] == "face");
tr2 = obj.problems(2).GlobalSpaces([obj.problems(2).GlobalSpaces.Support] == "face");
assert(isscalar(S2.S12) &&S2.d_g(1) == 3*d1, "InterfaceProblem:stokesBiotLayout", ...
    "stokes_biot needs one three-component trace on the Biot side " + ...
    "(LegendreBasis(k).to_vector(3)) and no other global space.");
n_loc2 = size(S2.S12{1}, 1);
assert(n_loc2 == dS + 6*d2, "InterfaceProblem:stokesBiotLayout", ...
    "Biot local layout [sigma, eta, gamma, q, p] expected (%d dofs), got %d.", ...
    dS + 6*d2, n_loc2);
G_fluid = get_affine_map(mesh_1, part_1.elements, inverse = true);

%% displacement history h* per patch, transported onto Sigma
h0 = c.history_skeleton;
if isscalar(h0)
    h0 = h0 + zeros(2*d1, n_p);
else
    h0 = h0(:, part_2.faces(:));
end
history_extension = c.history_extension;
if isempty(history_extension) && ~(isequal(c.history_stress, 0) ...
        && isequal(c.history_pressure, 0) && isequal(c.history_rotation, 0))
    history_extension = history_transport(c, k, mesh_2, part_2);
end
if ~isempty(history_extension)
    assert(isequal(size(history_extension), size(h0)), ...
        "InterfaceProblem:badHistoryExtension", ...
        "history extension must be %d x %d, got %s.", ...
        size(h0, 1), size(h0, 2), mat2str(size(history_extension)));
    h0 = h0 + history_extension;
end

%% interface sources at each side's own interface quadrature points
evaluate_on_points = @(f, x) cellfun(@(f) ...
    reshape(f(x(1, :, :), x(2, :, :)) + 0*x(1, :, :), sq, 1, 1, []), ...
    f, 'UniformOutput', 0);
T1 = evaluate_on_points(c.source_neumann, part_1.quadrature.interface);
T1 = cat(3, T1{:});
T2 = evaluate_on_points(c.source_dirichlet, part_2.quadrature.interface);
T2 = cat(3, T2{:});
T4 = evaluate_on_points(c.source_slip, part_2.quadrature.interface);
T4 = cat(3, T4{:});
x2 = part_2.quadrature.interface;
T3 = reshape(c.source_pressure(x2(1, :, :), x2(2, :, :)) + 0*x2(1, :, :), sq, []);

%% condensed blocks of every global space of each side
nt1 = numel(S1.S12);
S12f1 = cell(1, nt1);
for t = 1:nt1
    S12f1{t} = reshape(permute(S1.S12{t}, [1 2 7 5 3 4 6]), ...
        size(S1.S12{t}, 1), S1.d_g(t)*S1.f_g(t), []);
end
S12f2 = reshape(permute(S2.S12{1}, [1 2 7 5 3 4 6]), n_loc2, S2.d_g(1)*S2.f_g(1), []);

solid_sigma_basis = obj.problems(2).LocalSpaces(1).Basis;

[rows_acc, cols_acc, vals_acc] = deal({});
[beta_rows, beta_vals] = deal({});     % the level scalar's entries, kept apart
rhs_acc = zeros(n_dof, 1);
    function add(R, C, V)
        rows_acc{end+1} = repmat(R(:), 1, numel(C));
        cols_acc{end+1} = repmat(C(:).', numel(R), 1);
        vals_acc{end+1} = V;
    end

for e = 1:n_p
    K1 = part_1.elements(e);  F1 = part_1.faces(e);
    K2 = part_2.elements(e);  F2 = part_2.faces(e);
    fd1 = get_global_dof(tr1, subset=F1);  DOF1 = fd1(:);           % 2 d1
    fd2 = get_global_dof(tr2, subset=F2);  DOF2 = fd2(:) + off;     % 3 d1
    DOF2n = DOF2(1:d1);
    DOF2t = DOF2(d1 + (1:d1));
    DOF2u = DOF2(1:2*d1);
    DOF2p = DOF2(2*d1 + (1:d1));
    S12K1 = cell2mat(cellfun(@(A) A(:, :, K1), S12f1, 'UniformOutput', false));
    S12K2 = S12f2(:, :, K2);
    foot1 = cell2mat(cellfun(@(fp) fp(:, K1), S1.L.footprint(:).', 'UniformOutput', false).');
    foot2 = S2.L.footprint{1}(:, K2) + off;
    G_inner = G_fluid(:, :, e);

    nF = -reshape(part_1.quadrature.normals(:, :, e), 2, 1, []);
    nS = reshape(part_2.quadrature.normals(:, :, e), 2, 1, []);
    nS_hat = nS ./ vecnorm(nS, 2, 1);
    tS_hat = [-nS_hat(2, 1, :); nS_hat(1, 1, :)];
    nF_hat = nF ./ vecnorm(nF, 2, 1);

    %% polynomials
    s_fluid = part_1.quadrature.local(:, :, e);
    l1 = legendre_basis(k, s_fluid);
    L1 = vector_legendre_basis(k, s_fluid);
    s_solid = part_2.quadrature.local(:, :, e);
    l2 = legendre_basis(k, s_solid);
    L2 = vector_legendre_basis(k, s_solid);

    D1 = vector_dubiner_basis(k, part_1.quadrature.local_2d(1, :, e), ...
        part_1.quadrature.local_2d(2, :, e));
    D2 = vector_dubiner_basis(k, part_2.quadrature.local_2d(1, :, e), ...
        part_2.quadrature.local_2d(2, :, e));
    P2_face = dubiner_basis(k, part_2.quadrature.local_2d(1, :, e).', ...
        part_2.quadrature.local_2d(2, :, e).');                       % sq x d2
    P2_iface = dubiner_basis(k, part_2.quadrature.interface_2d(1, :, e).', ...
        part_2.quadrature.interface_2d(2, :, e).');                   % sq x d2

    % fluid stress / rotation along the paths and at Sigma
    s_conn_inner = reshape(part_1.quadrature.path_2d(1, :, e, :), sq, tq).';
    t_conn_inner = reshape(part_1.quadrature.path_2d(2, :, e, :), sq, tq).';
    TP_P1 = PEERS_basis(k, s_conn_inner(:), t_conn_inner(:), G_inner);
    TP_P1 = reshape(TP_P1, [tq sq size(TP_P1, 2:4)]);
    TP_A1 = asym_tensor_dubiner_basis(k, s_conn_inner(:), t_conn_inner(:));
    TP_A1 = reshape(TP_A1, [tq sq size(TP_A1, 2:4)]);

    s_iface_inner = part_1.quadrature.interface_2d(1, :, e);
    t_iface_inner = part_1.quadrature.interface_2d(2, :, e);
    E_P1 = PEERS_basis(k, s_iface_inner, t_iface_inner, G_inner);
    E_P1 = permute(E_P1, [3 4 1 2]);
    E_P1_N = permute(pagemtimes(E_P1, nF), [3 4 1 2]);
    E_P1_Nh = permute(pagemtimes(E_P1, nF_hat), [3 4 1 2]);
    E_A1 = asym_tensor_dubiner_basis(k, s_iface_inner, t_iface_inner);
    E_A1 = permute(E_A1, [3 4 1 2]);
    E_A1_N = permute(pagemtimes(E_A1, nF), [3 4 1 2]);

    % Biot stress / rotation / pressure along the paths and at Sigma
    s_conn_outer = reshape(part_2.quadrature.path_2d(1, :, e, :), sq, tq).';
    t_conn_outer = reshape(part_2.quadrature.path_2d(2, :, e, :), sq, tq).';
    TP_P2 = evaluate(solid_sigma_basis, [s_conn_outer(:).'; t_conn_outer(:).'], ...
        mesh_2, element_list = K2);
    TP_P2 = reshape(TP_P2, [tq sq size(TP_P2, 2:4)]);
    TP_A2 = asym_tensor_dubiner_basis(k, s_conn_outer(:), t_conn_outer(:));
    TP_A2 = reshape(TP_A2, [tq sq size(TP_A2, 2:4)]);
    TP_D2 = reshape(dubiner_basis(k, s_conn_outer(:), t_conn_outer(:)), tq, sq, d2);
    PI_p = reshape(sum(reshape(wTP, tq, 1, 1) .* TP_D2, 1), sq, d2);  % int_path phi

    s_iface_outer = part_2.quadrature.interface_2d(1, :, e);
    t_iface_outer = part_2.quadrature.interface_2d(2, :, e);
    E_P2 = evaluate(solid_sigma_basis, [s_iface_outer; t_iface_outer], ...
        mesh_2, element_list = K2);
    E_P2 = permute(E_P2, [3 4 1 2]);
    E_P2_N = permute(pagemtimes(E_P2, nS), [3 4 1 2]);
    E_A2 = asym_tensor_dubiner_basis(k, s_iface_outer, t_iface_outer);
    E_A2 = permute(E_A2, [3 4 1 2]);
    E_A2_N = permute(pagemtimes(E_A2, nS), [3 4 1 2]);

    % path integrals of the velocity / displacement reconstructions
    mF = reshape(part_1.quadrature.interface(:, :, e) ...
        - part_1.quadrature.physical(:, :, e), 2, 1, []);
    TP_P1 = permute(TP_P1, [4 5 2 1 3]);
    TP_P1_tr_Id = pagemtimes(TP_P1(1, 1, :, :, :) + TP_P1(2, 2, :, :, :), eye(2));
    TP_P1_grad_M = permute(pagemtimes(c1_1*TP_P1 - c2_1*TP_P1_tr_Id, mF), [3 4 5 1 2]);
    TP_INT_P1 = tensorprod(TP_P1_grad_M, wTP, 2, 1);
    TP_A1 = permute(TP_A1, [4 5 2 1 3]);
    TP_A1_M = permute(pagemtimes(TP_A1, mF), [3 4 5 1 2]);
    TP_INT_A1 = tensorprod(TP_A1_M, wTP, 2, 1);

    mS2 = part_2.quadrature.interface(:, :, e) - part_2.quadrature.physical(:, :, e);
    mS = reshape(mS2, 2, 1, []);
    TP_P2 = permute(TP_P2, [4 5 2 1 3]);
    TP_P2_tr_Id = pagemtimes(TP_P2(1, 1, :, :, :) + TP_P2(2, 2, :, :, :), eye(2));
    TP_P2_grad_M = permute(pagemtimes(c1_2*TP_P2 - c2_2*TP_P2_tr_Id, mS), [3 4 5 1 2]);
    TP_INT_P2 = tensorprod(TP_P2_grad_M, wTP, 2, 1);
    TP_A2 = permute(TP_A2, [4 5 2 1 3]);
    TP_A2_M = permute(pagemtimes(TP_A2, mS), [3 4 5 1 2]);
    TP_INT_A2 = tensorprod(TP_A2_M, wTP, 2, 1);
    % c_b p I in the displacement gradient: c_b int_path p m, [sq, d2, 2]
    TP_INT_p = c_b * cat(3, PI_p .* mS2(1, :).', PI_p .* mS2(2, :).');
    % Darcy pressure transfer: int_path (kappa^-1 q) . m on the q columns
    a_kappa = kappa.' \ mS2;                                        % 2 x sq
    PI_q = [PI_p .* a_kappa(1, :).', PI_p .* a_kappa(2, :).'];      % sq x 2 d2

    tau_1 = part_1.tau(e);
    tau_2 = part_2.tau(e);
    tau_p = tau_p_factor * tau_2;
    tau_s = c.stabilization_solid * tau_2;
    w_hat1 = tau_1*mu_f;
    Z_q = zeros(2*d1, 3*d2);     % the q, p columns of the Biot local block

    %% R_f: <Eh(sigma_f) + Eh(sigma_p), mu_1> = <T_1, mu_1>, fluid faces
    patch_len = part_1.area(e);
    integrate = @(D, L) patch_len * reshape(permute(pagemtimes(wS.*L, 'transpose', D, 'none'), [1 3 2]), 2*d1, []);
    S1_inner = [integrate(E_P1_N, l1), -tau_1*mu_f*integrate(D1, l1), ...
        integrate(0*E_A1_N, l1), integrate(0*D1(:, 1, :), l1)];
    S1_outer = [integrate(E_P2_N, l1), -tau_s*integrate(D2, l1), ...
        integrate(0*E_A2_N, l1), Z_q];
    add(DOF1, foot1, -S1_inner * S12K1);
    rhs_acc(DOF1) = rhs_acc(DOF1) - S1_inner * S1.S1f(:, :, 1, 1, K1);
    add(DOF1, foot2, -S1_outer * S12K2);
    rhs_acc(DOF1) = rhs_acc(DOF1) - S1_outer * S2.S1f(:, :, 1, 1, K2);
    add(DOF1, DOF1, integrate(L1, w_hat1 .* l1));
    add(DOF1, DOF2u, tau_s*integrate(L2, l1));
    rhs_acc(DOF1) = rhs_acc(DOF1) + integrate(T1(:, :, :, e), l1);

    %% rows tested on the Biot faces (scalar Legendre test, side-2 measure)
    patch_len = part_2.area(e);
    project = @(D, L, v) patch_len * ((wS(:).*L).' * ...
        (reshape(v(1, 1, :), [], 1).*reshape(D(:, :, 1), size(D, 1), []) ...
       + reshape(v(2, 1, :), [], 1).*reshape(D(:, :, 2), size(D, 1), [])));
    scalar = @(B) patch_len * ((wS(:).*l2).' * B);
    proj_n = @(D, L) project(D, L, nS_hat);
    proj_t = @(D, L) project(D, L, tS_hat);
    proj_f = @(D, L) project(D, L, nF_hat);
    Zp = zeros(d1, d2);

    %% R_n: <nu_p.(cdt Th(u_hat_f) - Th(eta_hat)) - cdt Eh(q_hat), l>
    Rn_inner = cdt*[proj_n(TP_INT_P1, l2), proj_n(0*D1, l2), ...
                    proj_n(TP_INT_A1, l2), proj_n(0*D1(:, 1, :), l2)];
    Eq_n = [P2_iface .* reshape(nS_hat(1, 1, :), [], 1), ...
            P2_iface .* reshape(nS_hat(2, 1, :), [], 1)];             % E(q).nu_p
    Rn_outer = [-[proj_n(TP_INT_P2, l2), proj_n(0*D2, l2), proj_n(TP_INT_A2, l2)], ...
        -cdt*scalar(Eq_n), -proj_n(TP_INT_p, l2) - cdt*tau_p*scalar(P2_face)];
    add(DOF2n, foot1, -Rn_inner * S12K1);
    rhs_acc(DOF2n) = rhs_acc(DOF2n) - Rn_inner * S1.S1f(:, :, 1, 1, K1);
    add(DOF2n, foot2, -Rn_outer * S12K2);
    rhs_acc(DOF2n) = rhs_acc(DOF2n) - Rn_outer * S2.S1f(:, :, 1, 1, K2);
    add(DOF2n, DOF1, cdt*proj_n(L1, l2));
    add(DOF2n, DOF2u, -proj_n(L2, l2));
    add(DOF2n, DOF2p, cdt*tau_p*scalar(l2));
    beta_rows{end+1} = DOF2n(:);  beta_vals{end+1} = -cdt*scalar(ones(sq, 1));
    rhs_acc(DOF2n) = rhs_acc(DOF2n) - proj_n(L2, l2)*h0(:, e) ...
        + cdt*proj_n(T2(:, :, :, e), l2);

    %% R_t: <t.Eh(sigma_f) + a t.(cdt Th(u_hat_f) - Th(eta_hat)), l>
    Rt_inner = [proj_t(E_P1_N, l2), -tau_1*mu_f*proj_t(D1, l2), ...
                proj_t(0*E_A1_N, l2), proj_t(0*D1(:, 1, :), l2)] ...
        + a_bjs*cdt*[proj_t(TP_INT_P1, l2), proj_t(0*D1, l2), ...
                     proj_t(TP_INT_A1, l2), proj_t(0*D1(:, 1, :), l2)];
    Rt_outer = -a_bjs*[proj_t(TP_INT_P2, l2), proj_t(0*D2, l2), proj_t(TP_INT_A2, l2), ...
        zeros(d1, 2*d2), proj_t(TP_INT_p, l2)];
    add(DOF2t, foot1, -Rt_inner * S12K1);
    rhs_acc(DOF2t) = rhs_acc(DOF2t) - Rt_inner * S1.S1f(:, :, 1, 1, K1);
    add(DOF2t, foot2, -Rt_outer * S12K2);
    rhs_acc(DOF2t) = rhs_acc(DOF2t) - Rt_outer * S2.S1f(:, :, 1, 1, K2);
    add(DOF2t, DOF1, proj_t(L1, w_hat1 .* l2) + a_bjs*cdt*proj_t(L1, l2));
    add(DOF2t, DOF2u, -a_bjs*proj_t(L2, l2));
    rhs_acc(DOF2t) = rhs_acc(DOF2t) - a_bjs*proj_t(L2, l2)*h0(:, e) ...
        + proj_t(T4(:, :, :, e), l2);

    %% R_p: <Tp(p_hat) + nu_f.Eh(sigma_f), l> = <T_3, l>
    Rp_inner = [proj_f(E_P1_Nh, l2), -tau_1*mu_f*proj_f(D1, l2), ...
                proj_f(0*E_A1_N, l2), proj_f(0*D1(:, 1, :), l2)];
    Rp_outer = [zeros(d1, dS + 3*d2), -scalar(PI_q), Zp];
    add(DOF2p, foot1, -Rp_inner * S12K1);
    rhs_acc(DOF2p) = rhs_acc(DOF2p) - Rp_inner * S1.S1f(:, :, 1, 1, K1);
    add(DOF2p, foot2, -Rp_outer * S12K2);
    rhs_acc(DOF2p) = rhs_acc(DOF2p) - Rp_outer * S2.S1f(:, :, 1, 1, K2);
    add(DOF2p, DOF1, w_hat1*proj_f(L1, l2));
    add(DOF2p, DOF2p, scalar(l2));
    rhs_acc(DOF2p) = rhs_acc(DOF2p) + scalar(T3(:, e));
end

%% pore-pressure level row: |Omega_P,h|^-1 (p, 1) on p = S1f - S12 x_hat
n_el2 = size(mesh_2.elements, 2);
[x_ref, w_ref] = get_2d_reference_quadrature(k + 1);
m_ref = DubinerBasis.evaluate_basis(k, x_ref(1, :), x_ref(2, :)).' * w_ref(:);   % d2 x 1
m_p = m_ref .* reshape(mesh_2.volume, 1, []) / sum(mesh_2.volume);           % d2 x n_el2
p_rows = dS + 5*d2 + (1:d2);
V = -pagemtimes(reshape(m_p, 1, d2, n_el2), S12f2(p_rows, :, :));            % 1 x n_foot x n_el2
foot_all = S2.L.footprint{1} + off;                                          % n_foot x n_el2
beta_rows{end+1} = a_dof*ones(numel(foot_all), 1);
beta_cols = reshape(foot_all, [], 1);
rhs_acc(a_dof) = rhs_acc(a_dof) + c.pore_pressure_mean ...
    - sum(m_p .* reshape(S2.S1f(p_rows, 1, 1, 1, :), d2, n_el2), "all");

rows = cellfun(@(v) v(:), rows_acc.', 'UniformOutput', false);
cols = cellfun(@(v) v(:), cols_acc.', 'UniformOutput', false);
vals = cellfun(@(v) v(:), vals_acc.', 'UniformOutput', false);
Hc = sparse(vertcat(rows{:}), vertcat(cols{:}), vertcat(vals{:}), n_dof, n_dof);
% the level scalar's column and row are assembled separately and added
% last, so the transmission entries are summed exactly as without it
n_col = sum(cellfun(@numel, beta_rows(1:end-1)));
Hc = Hc + sparse(vertcat(beta_rows{:}), ...
    [a_dof*ones(n_col, 1); beta_cols], ...
    [vertcat(beta_vals{:}); reshape(V, [], 1)], n_dof, n_dof);
Fc = rhs_acc;
end

%% ======================================================================
function ext = history_transport(c, k, mesh_2, part_2)
% Ext(eta*), (2*(k+1)) x n_patches, on the Biot side of the partition:
% trace_extension_coefficients of the effective stress sigma* + alpha p* I
% and the rotation gamma* (a scalar 0 is a zero field); the displacement
% gradient is the elastic compliance of the effective stress.
n_el = size(mesh_2.elements, 2);
d2 = nchoosek(k + 2, 2);
sigma_h = c.history_stress + zeros(PEERSBasis(k).Dimension, n_el);
p_h = c.history_pressure;
if ~isequal(p_h, 0)
    p_h = p_h + zeros(d2, n_el);
end
sigma_h = BiotProblem.effective_stress(sigma_h, p_h, c.biot_alpha, k);
gamma_h = c.history_rotation + zeros(d2, n_el);
solid = LinearElasticityProblem(k, young_modulus = c.young_modulus, ...
    poisson_ratio = c.poisson_ratio);
ext = trace_extension_coefficients(solid, mesh_2, part_2, sigma_h, gamma_h, Order = k);
end
