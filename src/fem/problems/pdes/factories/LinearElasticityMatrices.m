function [lhs, flux, rhs] = LinearElasticityMatrices(k, mesh, E, nu, options)
    % With density > 0 the inertial terms are added to the momentum
    % block, mirroring @LinearElasticityProblem/get_matrices:
    %   newmark: M/(beta dt^2) with the (u0 + dt v0)/(beta dt^2) and
    %            ((1/2-beta)/beta) a0 histories,
    %   onestep start (a state with u00 unset):
    %                      2M/dt^2 with 2(M u0)/dt^2 + 2(M v0)/dt,
    %   onestep marching:  M/dt^2 with M(2 u0 - u00)/dt^2.
    % The histories u0 = state.displacement, v0 = state.velocity,
    % a0 = state.acceleration, u00 = state.displacement_previous (the
    % solid sub-struct of @LinearElasticityProblem's state) are volume dof
    % arrays (2*d2 x n_el; [] or the scalar 0 = zero, but u00 = [] selects
    % the onestep start); rhs collects them as
    % a [-1, 2] block (empty when density = 0).
    arguments
        k double = 1
        mesh SimplicialGrid = UnitSquareMesh();
        E = 1.0;
        nu = 0.3333;
        options.density     (1,1) double = 0.0
        options.time_step   (1,1) double = 1.0
        options.time_scheme (1,1) string ...
            {mustBeMember(options.time_scheme, ["onestep", "newmark"])} = "onestep"
        options.newmark_beta (1,1) double = 1/4
        options.state                       struct = struct.empty
        % Winkler foundation b0 (u, v): the low-level path's hooke_modulus,
        % entering the momentum row's u-u block exactly as the inertial
        % mass does but without the time-scheme factor
        options.hooke_modulus (1,1) double = 0.0
        % Keep the element mean of tr(sigma) as a global unknown, as
        % IncompressibleElasticityMatrices does for the fluid: local ids
        % 1-4 = sigma, u, gamma, xi (DubinerBasis(0)), global ids 5-6 =
        % u_hat, rho (DubinerBasis(0) per element), with
        %   xi row:  (tr sigma, tr phi_1)_K - (phi_1, phi_1)_K rho_K = 0,
        %   rho row: <u_hat . n, phi_1>_dK - (phi_1, phi_1)_K rho_K / (2(lambda + mu)) = 0.
        % An exact reformulation: the local problems stay bounded as
        % lambda -> infinity and lambda enters only the 1/lambda rho
        % diagonal (no lambda-sized penalty in the condensed system).
        % false = the classical layout (local sigma, u, gamma; global u_hat).
        options.pressure_split (1,1) logical = false
    end
    lambda = E*nu/((1 + nu)*(1 - 2*nu));
    mu = E/2/(1 + nu);

    t = permute(mesh.tau, [3:6 2 1]);
    n = permute(mesh.exterior_unit_normals, [5 6 1 2 4 3]);

    S_k = PEERSBasis(k);
    P_k = DubinerBasis(k);
    L_k = LegendreBasis(k);
    
    c1 = 1 / (2 * mu);
    c2 = lambda / (4 * mu * (lambda + mu));
    
    M_S = peers_component_mass(S_k, mesh);
    d2 = P_k.Dimension;
    M_A = -M_S(:, d2 + (1:d2), 2, 1, :) + M_S(:, 2*d2 + (1:d2), 1, 2, :);

    M_T_S = pages(mass_matrix(mesh, S_k, S_k, ...
        trial_operator=TraceOperator(), test_operator=TraceOperator()));
    cv = convection_matrix(mesh, P_k, P_k);
    C = cat(3, pages(cv.x), pages(cv.y));
    PP = face_pages(boundary_matrix(mesh, P_k, P_k));
    LP = face_pages(boundary_matrix(mesh, P_k, L_k));
    n_LP = n .* LP;

    ll = pages(mass_matrix(mesh, L_k, L_k));
    LL = reshape(ll(:,:,:,:,mesh.boundaries.'), size(ll,1),size(ll,2),1,1,[],3);

    C_S = row_wise_operator(C);
    bubble = zeros(size(C_S, 1), S_k.Dimension - 4*d2, 1, 1, size(C_S, 5));
    C_S = [C_S, bubble];
    
    N_LS = row_wise_operator(n_LP);
    bubble = zeros(size(N_LS, 1), S_k.Dimension - 4*d2, 1, 1, size(N_LS, 5), size(N_LS, 6));
    N_LS = [N_LS, bubble];
    
    % (A sigma, tau) = c1 (sigma, tau) - c2 (tr sigma, tr tau):
    % sum the tensor component dims of the PEERS mass matrix
    sigma_tau_dx = c1*sum(M_S, 3:4);
    tr_sigma_tau_dx = c2*M_T_S;
    u_div_tau_dx = pagetranspose(C_S);
    gamma_tau_dx = M_A;
    % hybrid-trial blocks carry the face index in dim 7 (column faces),
    % as in DarcyMatrices, so the Schur complement couples all 3 faces
    U_t_n_dS = permute(pagetranspose(N_LS), [1:5 7 6]);
    split = options.pressure_split;
    iU = 4 + split;                                    % the u_hat id
    if split
        % phi_1 is the constant sigma_11 basis function, and
        % phi_1 + phi_(3 d2 + 1) = d0 I
        t_xi = M_T_S(:, 1, :, :, :);                   % (tr tau, tr phi_1)
        z_s = pagetranspose(t_xi);
        MS = sum(M_S, 3:4);
        z_rho = MS(1, 1, :, :, :);                     % (phi_1, phi_1)
        rho_row = sum(U_t_n_dS([1, 3*P_k.Dimension + 1], :, :, :, :, :, :), 1);   % <u_hat . n, d0>_dK
    end
    eq_1 = IndexedArray(sigma_tau_dx - tr_sigma_tau_dx, [1, 1]) ...
           + IndexedArray(u_div_tau_dx, [2, 1]) ...
           + IndexedArray(gamma_tau_dx, [3, 1]) ...
           + IndexedArray(-U_t_n_dS, [iU 1]);
    if split, eq_1 = eq_1 + IndexedArray(t_xi, [4 1]); end
    
    div_sigma_v_dx = -C_S;
    t_u_v_dS = sum(pagekron(eye(2), t .* PP), 6);
    if options.hooke_modulus ~= 0
        t_u_v_dS = t_u_v_dS + options.hooke_modulus ...
            * pagekron(eye(2), pages(mass_matrix(mesh, P_k, P_k)));
    end
    rhs = [];
    if options.density > 0
        n_el = size(mesh.elements, 2);
        mass = options.density * pagekron(eye(2), pages(mass_matrix(mesh, P_k, P_k)));
        dt = options.time_step;
        V_hist = zeros(2*d2, 1, 1, 1, n_el);
        % each history is a dof array or zero, read independently
        history = zeros(2*d2, n_el);
        h = solid_history(options, history);
        r = @(v) reshape(v, 2*d2, 1, 1, 1, []);
        if options.time_scheme == "newmark"
            beta = options.newmark_beta;
            t_u_v_dS = t_u_v_dS + mass/(beta*dt^2);
            u_0  = r(h.displacement);
            du_0 = r(h.velocity);
            V_hist = V_hist + pagemtimes(mass, u_0 + dt*du_0)/(beta*dt^2);
            a_0 = r(h.acceleration);
            V_hist = V_hist + (1/2 - beta)/beta * pagemtimes(mass, a_0);
        elseif ~isempty(options.state) && (~isfield(options.state, "displacement_previous") ...
                || isempty(options.state.displacement_previous))
            t_u_v_dS = t_u_v_dS + 2*mass/dt^2;
            u_0  = r(h.displacement);
            du_0 = r(h.velocity);
            V_hist = V_hist + 2*pagemtimes(mass, u_0)/dt^2 + 2*pagemtimes(mass, du_0)/dt;
        else
            t_u_v_dS = t_u_v_dS + mass/dt^2;
            u_00 = r(h.displacement_previous);
            u_0  = r(h.displacement);
            V_hist = V_hist + pagemtimes(mass, 2*u_0 - u_00)/dt^2;
        end
        rhs = IndexedMatrix(IndexedArray(V_hist, [-1, 2]));
    end
    gamma_v_dx = 0 * M_A(1:2*d2, :, :, :, :);
    t_U_v_dS = permute(pagekron(eye(2), t .* pagetranspose(LP)), [1:5 7 6]);
    eq_2 = IndexedArray(div_sigma_v_dx, [1, 2]) ...
           + IndexedArray(t_u_v_dS, [2, 2]) ...
           + IndexedArray(gamma_v_dx, [3, 2]) ... 
           + IndexedArray(-t_U_v_dS, [iU, 2]);
    
    sigma_eta_dx = pagetranspose(M_A);
    eq_3 = IndexedArray(sigma_eta_dx, [1 3]);
    
    lhs = eq_1 + eq_2 + eq_3;
    if split
        lhs = lhs + IndexedArray(z_s, [1 4]) + IndexedArray(-z_rho, [6 4]);
    end
    
    sigma_n_V_dx = N_LS;
    t_u_V_dx = pagekron(eye(2), t .* LP);
    t_U_V_dx = pagekron(eye(2), t .* LL);

    flux = IndexedArray(sigma_n_V_dx, [1, iU]) ...
           + IndexedArray(-t_u_V_dx, [2, iU]) ...
           + IndexedArray(t_U_V_dx, [iU, iU]);
    if split
        flux = flux + IndexedArray(rho_row, [iU, 6]) ...
                    + IndexedArray(-z_rho/(2*(lambda + mu)), [6, 6]);
    end
end

function h = solid_history(options, template)
    % The volume histories u0, v0, a0, u00 as arrays shaped like the
    % template, from options.state (empty state or missing field = zero).
    names = ["displacement", "velocity", "acceleration", "displacement_previous"];
    for i = 1:numel(names)
        value = [];
        if isfield(options.state, names(i)) && ~iscell(options.state.(names(i)))
            value = options.state.(names(i));
        end
        h.(names(i)) = default_zeros(value, template);
    end
end

function A = default_zeros(A, template)
    % Empty or scalar-0 histories expand to a zero array shaped like the template.
    if isempty(A) || isequal(A, 0.0)
        A = zeros(size(template));
    end
end

function M = peers_component_mass(basis, mesh)
% Component-paged PEERS Gram [d, d, 2, 2, n_el]: the (i, j) tensor
% component pages are kept uncontracted so M_A can slice the asymmetry
% blocks; sum(M, 3:4) recovers the full contraction (sigma, tau).
[x_local, w] = get_2d_reference_quadrature(basis.Order + 1);
vol = permute(mesh.volume, [1 3 4 5 2]);
p = basis.evaluate(x_local, mesh);
M = vol .* pagemtimes(w.*p, 'transpose', p, 'none');
end

function R_A = row_wise_operator(A)
    index_into = @(A, i) A(:, :, max(i, 1), :, :, :) .* permute(sign(i), [3 2 1]);
    resqueeze_cols = @(A) reshape(A, size(A, 1), [], 1, 1, size(A, 5), size(A, 6));
    resqueeze_rows = @(A) reshape(permute(A, [1 3 2 4 5 6]), [], size(A, 2), 1, 1, size(A, 5), size(A, 6));

    A_x = resqueeze_cols(index_into(A, [1; 0; 2; 0]));
    A_y = resqueeze_cols(index_into(A, [0; 1; 0; 2]));
    R_A = resqueeze_rows(cat(3, A_x, A_y));
end
