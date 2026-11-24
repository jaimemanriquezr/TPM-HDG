function sol = hdg_solve_cartesian(mesh, k, f, g, options)
%HDG_SOLVE_CARTESIAN  HDG for advection-diffusion with Q^k on a CartesianGrid.
%
%   sol = hdg_solve_cartesian(mesh, k, f, g, epsilon = 1, beta = [0; 0], ...
%                             tau = [], length_scale = 1, reaction = 0, ...
%                             mass_rhs = [])
%
%   Solves -epsilon lap(u) + beta . grad(u) + c u = f + w in the grid's
%   rectangle (c = reaction >= 0 constant; w = mass_rhs, a Q^k field given
%   by its Dimension-by-n_el QkBasis coefficients, or [] for w = 0),
%   u = g on its boundary, for constant epsilon > 0 and constant beta, in
%   the mixed form q + grad u = 0, div(epsilon q + beta u) + c u = f + w, with q in
%   [Q^k]^2, u in Q^k (QkBasis) on each element, the trace u_hat in P^k
%   (LegendreBasis) on each face, and the HDG flux
%       Fhat.n = epsilon q.n + (beta.n) u_hat + tau (u - u_hat).
%   Per element K, for all (r, w),
%       (q, r) - (u, div r) + <u_hat, r.n> = 0
%       -(epsilon q + beta u, grad w) + (c u, w) + <Fhat.n, w> = (f + w_h, w)
%   and on each interior face  sum_K <Fhat.n, mu> = 0. (q, u) are
%   condensed element by element; u_hat on boundary faces is the L2
%   projection of g.
%
%   One implicit Euler step of u_t - epsilon lap(u) + beta . grad(u) = f
%   is reaction = 1/dt, mass_rhs = u_old/dt (u_old = sol.u of the previous
%   step), f = f(t_new); both enter the local problems only.
%
%   tau is a scalar, or empty (default) for tau = epsilon/length_scale +
%   |beta.n| on each element side (single-valued on every face).
%   With epsilon = 1, beta = 0 and a scalar tau this is
%   hdg_poisson_cartesian.
%
%   f, g are function handles of (x, y). sol has fields
%     u      Dimension-by-n_el coefficients of u (QkBasis(k))
%     q      2*Dimension-by-n_el coefficients of q, [q_x; q_y]
%     u_hat  (k+1)-by-n_faces coefficients of u_hat (LegendreBasis(k))
%     n_trace_dof   number of global trace unknowns (interior faces)
%   and, with return_system = true, the condensed trace system:
%     K, G      sparse matrix and right-hand side over all face dof
%               (face f, mode m at (f-1)(k+1) + m), before the Dirichlet
%               elimination
%     K_local   4(k+1)-by-4(k+1)-by-n_el local condensed matrices, sides
%               S, E, N, W, modes fastest
%     load_local  4(k+1)-by-Dimension-by-n_el local load maps: the element's
%               contribution to G is load_local(:, :, K) * (f, phi_i)_K
%     qu_load_local, qu_trace_local   3 Dimension-by-Dimension-by-n_el and
%               3 Dimension-by-4(k+1)-by-n_el local recovery maps of
%               [q_x; q_y; u]_K = qu_load_local(:, :, K) * (f + w_h, phi_i)_K
%                     + qu_trace_local(:, :, K) * u_hat_K (sides S, E, N, W)
%     free, fixed   indices of the interior and Dirichlet dof
%
%   The local matrices are built from reference matrices scaled by hx, hy,
%   so they are exact for any (non-uniform) CartesianGrid.
%
%See also: hdg_poisson_cartesian, CartesianGrid, QkBasis, LegendreBasis
arguments
    mesh CartesianGrid
    k (1,1) double {mustBeNonnegative, mustBeInteger}
    f function_handle
    g function_handle
    options.epsilon (1,1) double {mustBePositive} = 1
    options.beta (2,1) double = [0; 0]
    options.tau double {mustBeScalarOrEmpty} = []
    options.length_scale (1,1) double {mustBePositive} = 1
    options.reaction (1,1) double {mustBeNonnegative} = 0
    options.mass_rhs double = []
    options.return_system (1,1) logical = false
end
epsilon = options.epsilon;
beta = options.beta;
basis = QkBasis(k);
d = basis.Dimension;            % (k+1)^2 element functions
m = k + 1;                      % face functions per face
n_el = size(mesh.elements, 2);
n_f = size(mesh.faces, 2);

% ---- reference matrices on [0, 1]^2 -----------------------------------
[x_ref, w_ref] = get_2d_reference_quadrature(k + 1, "Quad");
[P, Pr, Ps] = QkBasis.evaluate_basis(k, x_ref);
Dr = P.' * (w_ref .* Pr);       % Dr(a, b) = int phi_a d_r phi_b
Ds = P.' * (w_ref .* Ps);

[t, w_t] = get_1d_reference_quadrature(k + 1, "Quad");
mu = LegendreBasis.evaluate_basis(k, t.');      % n_t x m, orthonormal on [0, 1]
[E, F] = deal(cell(1, 4));
for side = 1:4
    Pside = QkBasis.evaluate_basis(k, mesh.get_2d_local_boundary(t(:), 1, side));
    E{side} = Pside.' * (w_t .* mu);            % int_side phi_i mu_m  (d x m)
    F{side} = Pside.' * (w_t .* Pside);         % int_side phi_i phi_j (d x d)
end
normal = [0, 1, 0, -1; -1, 0, 1, 0];            % outward normals of sides S, E, N, W
beta_n = beta.' * normal;                       % beta . n on sides S, E, N, W
if isempty(options.tau)
    tau = epsilon/options.length_scale + abs(beta_n);
else
    tau = options.tau * ones(1, 4);
end

% ---- local solvers -----------------------------------------------------
xq = mesh.get_physical_points(x_ref);           % 2 x n_q x n_el
f_q = reshape(f(xq(1, :, :), xq(2, :, :)), [], n_el);
F_rhs = mesh.volume .* (P.' * (w_ref .* f_q));  % (f, w_i), d x n_el
if ~isempty(options.mass_rhs)
    assert(isequal(size(options.mass_rhs), [d, n_el]), ...
        "hdg_solve_cartesian: mass_rhs must be Dimension-by-n_el")
    F_rhs = F_rhs + mesh.volume .* options.mass_rhs;   % (w_h, w_i), mass = |K| I
end

hx = mesh.hx; hy = mesh.hy;
[ix, jy] = ind2sub([mesh.Nx, mesh.Ny], 1:n_el);
n_loc = 4*m;
[K_loc, G_loc] = deal(zeros(n_loc, n_loc, n_el), zeros(n_loc, n_el));
[Ainv_B, Ainv_F] = deal(zeros(3*d, n_loc, n_el), zeros(3*d, n_el));
if options.return_system
    load_loc = zeros(n_loc, d, n_el);
    [qu_load_loc, qu_trace_loc] = deal(zeros(3*d, d, n_el), zeros(3*d, n_loc, n_el));
end
for K = 1:n_el
    h = [hx(ix(K)), hy(jy(K))];
    L = [h(1), h(2), h(1), h(2)];               % side lengths S, E, N, W
    vol = h(1)*h(2);

    Mq = vol*eye(2*d);
    Bq = [h(2)*Dr.'; h(1)*Ds.'];                % (u_j, div r_i), 2d x d
    G = beta(1)*h(2)*Dr.' + beta(2)*h(1)*Ds.';  % (phi_j, beta . grad phi_i)
    [Cq, Tuh, Buh] = deal(zeros(2*d, n_loc), zeros(d, n_loc), zeros(d, n_loc));
    Tuu = zeros(d);
    for side = 1:4
        cols = (side - 1)*m + (1:m);
        Cq(:, cols) = L(side) * [normal(1, side)*E{side}; normal(2, side)*E{side}];
        Tuh(:, cols) = tau(side)*L(side)*E{side};
        Buh(:, cols) = beta_n(side)*L(side)*E{side};
        Tuu = Tuu + tau(side)*L(side)*F{side};
    end
    Thh = kron(diag(tau.*L), eye(m));           % tau <mu, mu>, orthonormal mu
    Bhh = kron(diag(beta_n.*L), eye(m));        % (beta.n) <mu, mu>

    % local rows:  A (q; u) + B u_hat = (0; f)
    A = [Mq, -Bq; epsilon*Bq.', Tuu - G + options.reaction*vol*eye(d)];
    B = [Cq; Buh - Tuh];
    % transmission rows:  C'(q; u) + (Bhh - Thh) u_hat = 0
    C = [epsilon*Cq; Tuh];
    Ainv_B(:, :, K) = A \ B;
    Ainv_F(:, K) = A \ [zeros(2*d, 1); F_rhs(:, K)];
    K_loc(:, :, K) = Thh - Bhh + C.' * Ainv_B(:, :, K);
    G_loc(:, K) = C.' * Ainv_F(:, K);
    if options.return_system
        Ainv_W = A \ [zeros(2*d, d); eye(d)];
        load_loc(:, :, K) = C.' * Ainv_W;
        qu_load_loc(:, :, K) = Ainv_W;
        qu_trace_loc(:, :, K) = -Ainv_B(:, :, K);
    end
end

% ---- global trace system -------------------------------------------------
dof = reshape((mesh.boundaries(:).' - 1)*m + (1:m).', n_loc, n_el);   % local -> global
rows = repmat(dof, n_loc, 1);
cols = repelem(dof, n_loc, 1);
Kg = sparse(rows(:), cols(:), K_loc(:), m*n_f, m*n_f);
Gg = accumarray(dof(:), G_loc(:), [m*n_f, 1]);

% Dirichlet data: L2 projection of g on the boundary faces (orthonormal mu)
bdry = mesh.dirfaces;
v1 = mesh.coordinates(:, mesh.faces(1, bdry));
v2 = mesh.coordinates(:, mesh.faces(2, bdry));
xg = reshape(v1, 2, 1, []) + reshape(v2 - v1, 2, 1, []) .* t;    % 2 x n_t x n_b
g_t = reshape(g(xg(1, :, :), xg(2, :, :)), [], numel(bdry));
u_hat = zeros(m, n_f);
u_hat(:, bdry) = mu.' * (w_t .* g_t);

fixed = reshape((bdry - 1)*m + (1:m).', [], 1);
free = setdiff(1:m*n_f, fixed);
u_hat_vec = u_hat(:);
u_hat_vec(free) = Kg(free, free) \ (Gg(free) - Kg(free, fixed)*u_hat_vec(fixed));
u_hat = reshape(u_hat_vec, m, n_f);

% ---- local recovery ------------------------------------------------------
qu = Ainv_F - reshape(pagemtimes(Ainv_B, reshape(u_hat_vec(dof), n_loc, 1, n_el)), 3*d, n_el);
sol.q = qu(1:2*d, :);
sol.u = qu(2*d+1:end, :);
sol.u_hat = u_hat;
sol.n_trace_dof = numel(free);
if options.return_system
    sol.K = Kg;                     % condensed trace matrix, all faces
    sol.G = Gg;                     % its right-hand side
    sol.K_local = K_loc;            % per element, sides S, E, N, W
    sol.load_local = load_loc;      % (f, phi_i)_K -> element load
    sol.qu_load_local = qu_load_loc;    % (f + w_h, phi_i)_K -> (q; u)_K
    sol.qu_trace_local = qu_trace_loc;  % u_hat_K -> (q; u)_K
    sol.free = free;
    sol.fixed = fixed;
end
end
