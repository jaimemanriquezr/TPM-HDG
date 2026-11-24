function out = eulerian_hyperelastic_solid_march(k, N, amplitude, dt, n_steps, options)
%EULERIAN_HYPERELASTIC_SOLID_MARCH  Fully Eulerian hyperelastic solid on a moving snapped carve.
%
%   out = eulerian_hyperelastic_solid_march(k, N, A, dt, n_steps, rotation = 0, ...
%       law = "neo_hookean", fields = "standard", ...)
%
%   The solid alone, fully Eulerian, in the current configuration Omega(t):
%   the disc of the given radius and centre c(t) = c0 + A sin(2 pi t) e_x.
%   At every t_n the solid mesh is the inside of the circle Sigma(t_n)
%   carved from the background StructuredSimplicialGrid(N, N, "crossed")
%   on [0, 1]^2 (cut elements discarded; only the cells near the disc are
%   built) with its interface vertices snapped onto Sigma(t_n)
%   (snap_to_curve, normal_sign = -1); every boundary face is Dirichlet.
%   Unknowns: the velocity v and the Eulerian displacement u, x - u(x, t)
%   = X the reference map, with
%     d_t u + (v.grad) u - v = g,   rho (d_t v + (v.grad) v) - div sigma = f,
%     sigma = sigma_E(grad u),
%   rho constant and g = 0 physically (a manufactured source otherwise).
%   law = "neo_hookean": with a = I - G = F^-1, j = det a, B = (a' a)^-1,
%     sigma_E(G) = j (mu (B - I) - lambda log(j) I),
%   the Cauchy stress of P = mu (F - F^-T) + lambda log J F^-T; law =
%   "svk" (Saint Venant-Kirchhoff): sigma_E(G) = j F S F', F = a^-1, S =
%   lambda tr(E) I + 2 mu E, E = (F' F - I)/2; law = "linear": sigma_E(G)
%   = C eps(G). Split law sigma_E(G) = C eps(G) +
%   N_E(G), N_E quadratic at 0 and symmetric.
%
%   Discretization (degree k): per element sigma in the PEERS space S_h,
%   u, v in P_k^2, gamma in P_k (skew multiplier); per face u_hat, v_hat in
%   P_k(F)^2. L_h is the PEERS lifting of (u, u_hat),
%   (L_h, tau)_K = -(u, div tau)_K + <u_hat, tau n>_dK. One BDF2 step on
%   the carve of t_n+1, c = 3/(2 dt), Z = (4 X~^n - X~^n-1)/3, beta the
%   transport velocity in K, beta_hat on the faces, n the outward normal:
%    (C)  (A sigma, tau) + (u, div tau) + (gamma, tau) - <u_hat, tau n> - (A N_E(L_h), tau) = 0
%    (S)  (sigma, eta) - (N_E(L_h), eta) = 0
%    (M)  c rho (v - Z_v, w) + rho (beta.grad v, w) - (div sigma, w)
%           + <(tau_v + rho (a_v - beta_hat.n)) (v - v_hat), w> = (f, w)
%    (K)  c (u - Z_u, psi) + (beta.grad u, psi) + <(a_u - beta_hat.n) (u - u_hat), psi>
%           - (v, psi) = (g, psi)
%    (Fv) sum_K <sigma n - (tau_v + rho a_v) (v - v_hat), mu>_dK = 0   (interior faces)
%    (Fu) sum_K <a_u (u - u_hat), mu>_dK = 0                           (interior faces)
%   with a_v = |beta_hat.n| (upwind flux of v), a_u = (beta_hat.n +
%   sqrt((beta_hat.n)^2 + delta^2))/2 + tau_u (smoothed upwind trace of u
%   plus displacement stabilization); defaults tau_v = rho c_s, tau_u =
%   c_s, delta = upwind_smoothing = 0.1 c_s, c_s = sqrt(E/rho). The
%   convection is the advective form, which annihilates constants. (C) and
%   (S) are the constitutive and weak-symmetry rows of HyperelasticProblem
%   (material "neo_hookean_eulerian" or "saint_venant_kirchhoff_eulerian";
%   for the linear law those of
%   LinearElasticityProblem). Boundary faces: v_hat = P g_v (face L2
%   projection of the exact velocity on the chord); displacement_boundary
%   = "dirichlet": u_hat = P g_u, "free": <u - u_hat, mu>_F = 0 (u_hat the
%   trace of u), "characteristic": <a_u (u - u_hat) + a_in (g_u - u_hat),
%   mu>_F = 0 with a_in = (-beta_hat.n + sqrt((beta_hat.n)^2 + delta^2))/2
%   (u_hat the outgoing trace on outflow, the data g_u enter on inflow
%   only; a coupled march may pass data.kinematic = struct(Zu, dt) to the
%   local operators, and g_u is then the BDF2 step of the kinematic
%   equation at the face points, Z_u + (2 dt/3)(v - (grad u) v), from the
%   iterate, differentiated in Newton). The face terms use one Gauss rule exact for degree 3k + 1
%   in the element rows and the trace rows.
%
%   Nonlinear iteration, one loop: each iteration is one linear solve for
%   the new iterate, with N_E linearized at L_h(u0, u_hat0) (Newton, the
%   class's get_matrices with Data.iterate) and beta = v0, beta_hat =
%   v_hat0 (Dirichlet data on boundary faces), a_v, a_u from beta_hat;
%   nonlinear = "newton" adds the derivatives of the transport terms with
%   respect to beta, beta_hat and the weights (d|b|/db = sign(b), the
%   smoothed a_u and a_in exactly), so that the loop is Newton's method;
%   "picard" lags them. Initial iterate: volume dofs 2 X~^n - X~^n-1 (also at
%   the first step, from the two initial levels), skeletons the face L2
%   projection of the mean of the adjacent element traces (on boundary
%   faces the inner trace), then the data on the Dirichlet v_hat dofs
%   (u_hat keeps the inner trace: the data of t_n+1 against the
%   extrapolated u would put a jump of order dt^2 |d_tt u| k^2/h into
%   L_h, enough to make det(I - L_h) negative under rotation at k = 3).
%   Stops when the
%   relative change of all dofs (v, u, v_hat, u_hat) is at most
%   picard_tolerance, at most picard_max_iterations iterations.
%
%   Moving carve: the time derivative is at fixed x, so a level X^m =
%   (v^m, u^m) on carve(t_m) is carried to carve(t_n+1) by the exact
%   two-mesh L2 projection with zero shift
%   (LinearElasticityProblem.translated_field_between), both BDF2 levels
%   directly from their own carves; on the part of a new element the old
%   carve does not cover, the degree-extension_degree Taylor polynomial of
%   the source old element's polynomial ([], the default: the full
%   polynomial).
%
%   Exact fields, rigid motion c(t), theta(t) = rotation sin(2 pi t), Q =
%   Q(theta):
%   fields = "standard": X(x, t) = Y - w(Y, t), Y = c0 + Q' (x - c(t)),
%     u = x - X, v = (I - grad u)^-1 d_t u, g = 0 (w = 0: the rigid motion);
%   fields = "quadratic": u, v quadratic in t, smooth in x, unrelated;
%   fields = "patch": u, v polynomials of total degree patch_degree in x
%     times (1 + t^2) (u) and (1 + t) (v); the rotation does not enter.
%   sigma = sigma_E(grad u), f and g from the equations. Start: the L2
%   projections of the exact (v, u) at t = 0 on carve(0) and at t = -dt
%   on carve(-dt).
%   initial_perturbation = [seed, amplitude]: a second state through the
%   same carves, the first state's initial levels plus amplitude times a
%   random (dv, du) of unit linearized energy on carve(0) (its level -1
%   copy that perturbation projected to carve(-dt)).
%
%   out: e_v, e_u, e_sigma (L2 errors at T, sigma the PEERS field against
%   the exact Cauchy stress), n_v, n_u, n_sigma (exact norms); step_errors
%   (n_steps x 3, [v, u, sigma], record_errors only); iterations and
%   converged (n_steps x 1; converged covers both states); uncovered
%   (n_steps x 2: the largest fraction of a new element and the total area
%   not covered by the old carve of t_n); deviation_energy (n_steps + 1):
%   1/2 rho |v_h - v|^2 + 1/2 (C eps(L_h - grad u), eps(L_h - grad u));
%   perturbation_energy (n_steps + 1): the linearized energy
%   1/2 rho |dv|^2 + 1/2 (C eps(dL), eps(dL)) of the difference of the two
%   states (dL the lifting of (du, du_hat)); field_bounds = [max |grad u|
%   (Frobenius), max |v|/c_s, min det(I - grad u)] of the exact fields
%   over the run (quadrature points of every carve); the final v, u, sigma
%   (element dofs), u_hat, v_hat (2*(k+1) x n_faces), mesh; seconds;
%   changes (per step, the relative change of every iteration);
%   jacobian_error (jacobian_check = true): at the first step's initial
%   iterate, the relative difference between the central difference of the
%   residual in a random direction and the Newton operator applied to it,
%   for steps 1e-3, ..., 1e-8 (relative to the iterate's norm).
%   operators = true: no march; out holds prm, hp, the exact fields, the
%   carve and handles to the local functions a coupled march needs
%   (eulerian_hyperelastic_fsi_march).
%
%See also: velocity_stress_solid_march, HyperelasticProblem,
%   LinearElasticityProblem.translated_field_between
arguments
    k (1,1) double
    N (1,1) double
    amplitude (1,1) double
    dt (1,1) double
    n_steps (1,1) double
    options.rotation (1,1) double = 0
    options.law (1,1) string {mustBeMember(options.law, ["neo_hookean", "linear", "svk"])} = "neo_hookean"
    options.fields (1,1) string {mustBeMember(options.fields, ["standard", "quadratic", "patch"])} = "standard"
    options.patch_degree double = []
    options.w (2,1) string = ["0.02*sin(5*X + 3*Y)*cos(3*t)"; "0.02*cos(4*X - 5*Y)*(1 + t^2)"]
    options.displacement_boundary (1,1) string {mustBeMember(options.displacement_boundary, ...
        ["dirichlet", "free", "characteristic"])} = "dirichlet"
    options.extension_degree double = []
    options.tau_v double = []
    options.tau_u double = []
    options.upwind_smoothing double = []
    options.nonlinear (1,1) string {mustBeMember(options.nonlinear, ["newton", "picard"])} = "newton"
    options.picard_tolerance (1,1) double = 1e-10
    options.picard_max_iterations (1,1) double = 30
    options.jacobian_check (1,1) logical = false
    options.operators (1,1) logical = false
    options.initial_perturbation double = []
    options.record_errors (1,1) logical = false
    options.center (2,1) double = [0.5; 0.5]
    options.radius (1,1) double = 0.2
    options.snap (1,1) logical = true
    options.rho (1,1) double = 1
    options.E (1,1) double = 10
    options.nu (1,1) double = 0.3
end
wall = tic;
rho = options.rho;
c_s = sqrt(options.E/rho);
prm.rho = rho;
prm.tau_v = default_to(options.tau_v, rho*c_s);
prm.tau_u = default_to(options.tau_u, c_s);
prm.delta = default_to(options.upwind_smoothing, 0.1*c_s);
prm.free = ismember(options.displacement_boundary, ["free", "characteristic"]);
prm.char = options.displacement_boundary == "characteristic";
prm.newton = options.nonlinear == "newton";
prm.c = 3/(2*dt);
[lambda, mu] = LinearElasticityProblem.get_lame_parameters(options.E, options.nu);
prm.lambda = lambda;
prm.mu = mu;
if isempty(options.patch_degree), options.patch_degree = k; end
cen = @(t) options.center + [amplitude*sin(2*pi*t); 0];
moving = amplitude ~= 0;   % a rotation does not change the carve
carve = @(t) carve_solid(local_background(N, cen(t), options.radius), ...
    CircularCurve(cen(t), options.radius), options.snap);
ex = exact_fields(options, amplitude, lambda, mu);
hp = [];
if options.law ~= "linear"
    material = struct("neo_hookean", "neo_hookean_eulerian", "svk", "saint_venant_kirchhoff_eulerian");
    hp = HyperelasticProblem(k, material = material.(options.law), young_modulus = options.E, ...
        poisson_ratio = options.nu);
end
sp0 = LinearElasticityProblem(k);   % translated_field_between and the boundary data
project_between = @(old, new, U) project_levels(sp0, old, new, k, U, options.extension_degree);
if options.operators
    % the pieces a coupled march builds its solid from (n_steps unused)
    out = struct("prm", prm, "hp", hp, "exact", ex, "carve", carve, ...
        "geometry", @(mesh) geometry(mesh, k, options), "local_operators", @local_operators, ...
        "project_between", project_between, "project_vector", @project_vector, ...
        "moments", @moments, "initial_skeleton", @initial_skeleton_from, "errors", @errors, ...
        "deviation", @deviation, "linearized_energy", @linearized_energy, ...
        "jacobian_check", @jacobian_check);
    return
end

% initial levels: exact (v, u) at t = 0 on carve(0), at t = -dt on carve(-dt)
mesh = carve(0);
geo = geometry(mesh, k, options);
mesh_prev = mesh;
if moving, mesh_prev = carve(-dt); end
level = project_vector(mesh, k, @(x, y) [ex.v(x, y, 0); ex.u(x, y, 0)]);
level_prev = project_vector(mesh_prev, k, @(x, y) [ex.v(x, y, -dt); ex.u(x, y, -dt)]);
lam = initial_skeleton(geo, k, level, ex, 0, prm);
twin = ~isempty(options.initial_perturbation);
if twin
    stream = RandStream("mt19937ar", Seed = options.initial_perturbation(1));
    d = randn(stream, size(level));
    d_lam = initial_skeleton(geo, k, d, [], 0, prm);
    d = d/sqrt(linearized_energy(geo, k, prm, d, d_lam));
    d = options.initial_perturbation(2)*d;
    level2 = level + d;
    if moving
        level2_prev = level_prev + project_between(mesh, mesh_prev, d);
    else
        level2_prev = level_prev + d;
    end
    lam2 = initial_skeleton(geo, k, level2, ex, 0, prm);
end

deviation_energy = nan(n_steps + 1, 1);
perturbation_energy = nan(n_steps + 1, 1);
deviation_energy(1) = deviation(geo, k, prm, level, lam, ex, 0);
if twin, perturbation_energy(1) = linearized_energy(geo, k, prm, level2 - level, lam2 - lam); end
bounds = field_bounds(mesh, k, ex, 0, c_s);
iterations = zeros(n_steps, 1);
converged = true(n_steps, 1);
uncovered = zeros(n_steps, 2);
step_errors = nan(n_steps, 3);
changes = cell(n_steps, 1);
jacobian_error = [];
for n = 1:n_steps
    t1 = n*dt;
    if moving
        mesh_new = carve(t1);
        geo_new = geometry(mesh_new, k, options);
        [Xn, info] = project_between(mesh, mesh_new, level);
        Xm = project_between(mesh_prev, mesh_new, level_prev);
        uncovered(n, :) = [max(info.uncovered), sum(info.uncovered(:).*mesh_new.volume(:))];
        if twin
            Xn2 = project_between(mesh, mesh_new, level2);
            Xm2 = project_between(mesh_prev, mesh_new, level2_prev);
        end
        mesh_prev = mesh;
        mesh = mesh_new;
        geo = geo_new;
    else
        [Xn, Xm] = deal(level, level_prev);
        if twin, [Xn2, Xm2] = deal(level2, level2_prev); end
    end
    data = step_data(geo, k, ex, t1, prm);
    if options.jacobian_check && n == 1
        level0 = 2*Xn - Xm;
        jacobian_error = jacobian_check(geo, k, hp, prm, data, (4*Xn - Xm)/3, level0, ...
            initial_skeleton_from(geo, k, level0, data));
    end
    [level_new, lam, S, iterations(n), converged(n), changes{n}] = bdf2_step(geo, k, hp, prm, data, ...
        Xn, Xm, n, options);
    [level_prev, level] = deal(level, level_new);
    if twin
        [level2_new, lam2, ~, ~, ok] = bdf2_step(geo, k, hp, prm, data, Xn2, Xm2, n, options);
        converged(n) = converged(n) && ok;
        [level2_prev, level2] = deal(level2, level2_new);
        perturbation_energy(n + 1) = linearized_energy(geo, k, prm, level2 - level, lam2 - lam);
    end
    deviation_energy(n + 1) = deviation(geo, k, prm, level, lam, ex, t1);
    b = field_bounds(mesh, k, ex, t1, c_s);
    bounds = [max(bounds(1:2), b(1:2)), min(bounds(3), b(3))];
    if options.record_errors
        [e, ~] = errors(geo, k, level, S, ex, t1);
        step_errors(n, :) = e;
    end
end
T = n_steps*dt;
[e, nrm] = errors(geo, k, level, S, ex, T);
d2 = nchoosek(k + 2, 2);
n_dof = (k + 1)*size(mesh.faces, 2);
hat = @(a) [reshape(a(1:n_dof), k + 1, []); reshape(a(n_dof + 1:end), k + 1, [])];
out = struct("e_v", e(1), "e_u", e(2), "e_sigma", e(3), "n_v", nrm(1), "n_u", nrm(2), ...
    "n_sigma", nrm(3), "step_errors", step_errors, "iterations", iterations, ...
    "converged", converged, "uncovered", uncovered, "deviation_energy", deviation_energy, ...
    "perturbation_energy", perturbation_energy, "field_bounds", bounds, ...
    "v", level(1:2*d2, :), "u", level(2*d2 + 1:end, :), "sigma", S, ...
    "u_hat", hat(lam(1:2*n_dof)), "v_hat", hat(lam(2*n_dof + 1:end)), ...
    "mesh", mesh, "n_elements", size(mesh.elements, 2), "seconds", toc(wall));
out.changes = changes;
out.jacobian_error = jacobian_error;
end

function x = default_to(x, value)
if isempty(x), x = value; end
end

function grid = local_background(N, center, radius)
% the cells of StructuredSimplicialGrid(N, N, "crossed") on [0, 1]^2 that
% lie within two cells of the disc's bounding box, as their own crossed
% grid (the same triangles; the carve reads only those near the disc)
lo = max(floor((center - radius)*N) - 2, 0);
hi = min(ceil((center + radius)*N) + 2, N);
grid = StructuredSimplicialGrid(hi(1) - lo(1), hi(2) - lo(2), lo(1)/N, hi(1)/N, lo(2)/N, hi(2)/N, ...
    style = "crossed");
end

function mesh = carve_solid(background, curve, snap)
% the inside of the circle, cut elements discarded, interface vertices
% snapped onto it; every boundary face Dirichlet; tau = 1 (the class rows
% taken from get_matrices do not depend on it)
if snap
    [~, mesh] = background.split(curve, max_normal_deviation = 90);
    mesh = snap_to_curve(mesh, curve, normal_sign = -1);
else
    [~, mesh] = background.split(curve);
end
is_boundary = accumarray(mesh.boundaries(:), 1, [size(mesh.faces, 2), 1]) == 1;
mesh.faces(3, is_boundary) = mesh.boundary_tags("Dirichlet");
mesh = add_stabilization(mesh, 1);
end

function [V, info] = project_levels(sp0, old, new, k, U, q)
% the levels [v; u] of old on new, zero shift, degree-q extension
if isempty(q)
    [V, info] = sp0.translated_field_between(old, new, k, k, U, [0; 0]);
else
    [V, info] = sp0.translated_field_between(old, new, k, k, U, [0; 0], extension_degree = q);
end
end

function ex = exact_fields(options, amplitude, lambda, mu)
% handles of (x, y, t): v, u, f, g (2 x n), sigma, grad u (rows 11, 21, 12, 22)
syms x y t real
xx = [x; y];
c = options.center + [amplitude*sin(2*pi*t); 0];
th = options.rotation*sin(2*pi*t);
Q = [cos(th), -sin(th); sin(th), cos(th)];
switch options.fields
    case "standard"
        syms X Y real
        Yr = options.center + Q.'*(xx - c);
        u = xx - (Yr - subs(str2sym(options.w(:)), [X; Y], Yr));
        G = jacobian(u, xx);
        v = (eye(2) - G) \ diff(u, t);
        g = sym([0; 0]);
    case "quadratic"
        u = 0.03*[sin(2*x + y); cos(x - 2*y)] + t*0.02*[cos(3*y); sin(2*x)] ...
            + t^2*0.02*[x*y; x - y^2];
        v = 0.2*[cos(x + 2*y); sin(2*x - y)] + t*0.1*[sin(y); cos(x)] + t^2*0.1*[x*y; x^2];
        G = jacobian(u, xx);
        g = diff(u, t) + G*v - v;
    case "patch"
        U = sym([0; 0]);
        V = sym([0; 0]);
        i = 0;
        for a = 0:options.patch_degree
            for b = 0:options.patch_degree - a
                i = i + 1;
                m = x^a*y^b;
                V = V + m*[cos(1 + i); sin(2 + 3*i)];
                U = U + m*[cos(3 + i); sin(5 + 2*i)];
            end
        end
        u = 0.05*U*(1 + t^2);
        v = V*(1 + t);
        G = jacobian(u, xx);
        g = diff(u, t) + G*v - v;
end
if options.law == "linear"
    sigma = mu*(G + G.') + lambda*trace(G)*eye(2);
elseif options.law == "svk"
    a = eye(2) - G;
    Fd = inv(a);
    Eg = (Fd.'*Fd - eye(2))/2;
    sigma = det(a)*Fd*(lambda*trace(Eg)*eye(2) + 2*mu*Eg)*Fd.';
else
    a = eye(2) - G;
    j = det(a);
    sigma = j*(mu*(inv(a.'*a) - eye(2)) - lambda*log(j)*eye(2));
end
f = options.rho*(diff(v, t) + jacobian(v, xx)*v) ...
    - [diff(sigma(1, 1), x) + diff(sigma(1, 2), y); diff(sigma(2, 1), x) + diff(sigma(2, 2), y)];
ex.v = stacked(v, x, y, t);
ex.u = stacked(u, x, y, t);
ex.f = stacked(f, x, y, t);
ex.g = stacked(g, x, y, t);
ex.sigma = stacked(sigma(:), x, y, t);
ex.grad_u = stacked(G(:), x, y, t);
end

function h = stacked(e, x, y, t)
% handle (x, y, t) -> [e_1(x, y, t); e_2; ...] for column vectors x, y
% of equal size (rows of the result: component blocks)
hs = arrayfun(@(ei) matlabFunction(ei, Vars = [x, y, t], Optimize = false), e, UniformOutput = false);
h = @(xx, yy, tt) cell2mat(cellfun(@(f) f(xx, yy, tt) + 0*xx, hs, UniformOutput = false));
end

function [rs, w] = rule(q)
[rs, w] = get_2d_reference_quadrature(q);
w = w(:)/sum(w);
end

function [xs, ys] = points(mesh, rs)
% physical points (nq x n_el each) of the reference points rs
n_el = size(mesh.elements, 2);
xy = mesh.get_2d_physical(repmat(rs, 1, 1, n_el), 1:n_el);
xs = reshape(xy(1, :, :), size(rs, 2), []);
ys = reshape(xy(2, :, :), size(rs, 2), []);
end

function V = project_vector(mesh, k, h)
% L2 projection onto P_k of every component block of h(x, y)
[rs, w] = rule(k + 3);
[xs, ys] = points(mesh, rs);
nq = numel(w);
s = h(xs(:), ys(:));
D = dubiner_basis(k, rs);
fit = (D.'*(w.*D)) \ (D.'.*w.');
n_c = size(s, 1)/numel(xs);
V = zeros(n_c*size(D, 2), size(xs, 2));
for i = 1:n_c
    V((i - 1)*size(D, 2) + (1:size(D, 2)), :) = fit*reshape(s((i - 1)*numel(xs) + (1:numel(xs))), nq, []);
end
end

function F = moments(mesh, k, h)
% (h_c, P_i)_K of a vector field h(x, y) (2 components), rule exact for
% degree 3k + 3: 2*d2 x 1 x n_el
[rs, w] = rule(ceil(3*k/2) + 1);
[xs, ys] = points(mesh, rs);
s = h(xs(:), ys(:));
D = dubiner_basis(k, rs);
n = numel(xs);
F = [(w.*D).'*reshape(s(1:n), numel(w), []); (w.*D).'*reshape(s(n + 1:2*n), numel(w), [])];
F = reshape(F.*mesh.volume(:).', [], 1, size(xs, 2));
end

function geo = geometry(mesh, k, options)
% the operators of one carve that do not depend on the iterate
geo.mesh = mesh;
n_el = size(mesh.elements, 2);
d1 = k + 1;
d2 = nchoosek(k + 2, 2);
lin = get_matrices(LinearElasticityProblem(k, young_modulus = options.E, ...
    poisson_ratio = options.nu), mesh, k);
geo.lin = lin;
d2peers = size(lin{1}, 1) - 3*d2;
geo.d2peers = d2peers;
geo.tU = lin{1}(1:d2peers, d2peers + (1:2*d2), :);   % (u, div tau)
geo.M = mass_matrix(mesh, k);
[~, nxLP, nyLP] = boundary_matrix(mesh, k);
O = zeros(3*d1, d2, n_el);
o = zeros(3*d1, d2peers - 4*d2, n_el);
geo.Sn = [nxLP, O, nyLP, O, o; O, nxLP, O, nyLP, o];   % <sigma n, mu>
[geo.R_u, geo.R_hat] = HyperelasticProblem.discrete_gradient_operators(mesh, k);
% face rule (Gauss, exact for degree 3k + 1) at element-local nodes:
% volume basis Pf, face basis Lf (permutation of the face parametrization
% as in convective_matrix), weights wa (with the face length), unit
% outward normals, boundary flags
[s, w1] = get_1d_reference_quadrature(ceil(3*k/2));
nq = numel(w1);
q2 = [[0*s(:), 1 - s(:)]; [1 - s(:), 0*s(:)]; [1 - s(:), s(:)]];
geo.Pf = dubiner_basis(k, q2(:, 1), q2(:, 2));
geo.q2 = q2;
% physical gradients of the volume basis at the face points (3nq x d2 x n_el)
[~, Pr, Ps] = dubiner_basis(k, q2(:, 1).', q2(:, 2).');
F = get_affine_map(mesh, output = "array");
dF = F(1, 1, :).*F(2, 2, :) - F(1, 2, :).*F(2, 1, :);
geo.Pfx = (F(2, 2, :).*Pr - F(2, 1, :).*Ps)./dF;
geo.Pfy = (-F(1, 2, :).*Pr + F(1, 1, :).*Ps)./dF;
l = legendre_basis(k, [s(:); 1 - s(:)]);
geo.Lf = zeros(3*nq, 3*d1, n_el);
area = mesh.area(mesh.boundaries);
on = mesh.normals.*kron(mesh.orientation, [1; 1]);
is_boundary = accumarray(mesh.boundaries(:), 1, [size(mesh.faces, 2), 1]) == 1;
[geo.wa, geo.nx, geo.ny, geo.bnd] = deal(zeros(3*nq, n_el));
for e = 1:3
    nodes = (e - 1)*nq + (1:nq);
    for p = 1:2
        sel = mesh.permutation(e, :) == p;
        geo.Lf(nodes, (e - 1)*d1 + (1:d1), sel) = repmat(l((p - 1)*nq + (1:nq), :), 1, 1, nnz(sel));
    end
    geo.wa(nodes, :) = w1(:).*area(e, :);
    geo.nx(nodes, :) = repmat(on(2*e - 1, :)./area(e, :), nq, 1);
    geo.ny(nodes, :) = repmat(on(2*e, :)./area(e, :), nq, 1);
    geo.bnd(nodes, :) = repmat(is_boundary(mesh.boundaries(e, :)).', nq, 1);
end
% global indices of the element hybrid dofs [u_hat x; u_hat y; v_hat x; v_hat y]
n_dof = d1*size(mesh.faces, 2);
face_dof = reshape(((mesh.boundaries(:) - 1)*d1 + (1:d1)).', 3*d1, n_el);
geo.idx = [face_dof; face_dof + n_dof; face_dof + 2*n_dof; face_dof + 3*n_dof];
geo.n_dof = n_dof;
end

function B = blk2(A)
% [A 0; 0 A] page-wise
Z = zeros(size(A, 1), size(A, 2), size(A, 3));
B = [A, Z; Z, A];
end

function [PP, LP, LL] = weighted_face(geo, omega)
% <omega P_j, P_i>, <omega P_j, L_i>, <omega L_j, L_i> over dK
n_el = size(omega, 2);
W = reshape(geo.wa.*omega, [], 1, n_el);
WP = W.*geo.Pf;
PP = pagemtimes(geo.Pf, "transpose", WP, "none");
LP = pagemtimes(geo.Lf, "transpose", WP, "none");
LL = pagemtimes(geo.Lf, "transpose", W.*geo.Lf, "none");
end

function data = step_data(geo, k, ex, t, prm)
% sources and boundary data at t on the carve of geo
mesh = geo.mesh;
data.f = moments(mesh, k, @(x, y) ex.f(x, y, t));
data.g = moments(mesh, k, @(x, y) ex.g(x, y, t));
[data.lam, data.dir] = boundary_values(geo, k, ex, t, prm);
if prm.char
    data.gu_face = face_values(geo, @(x, y) ex.u(x, y, t));
end
end

function V = face_values(geo, h)
% values of a 2-component field h(x, y) at the face quadrature points of
% every element: 3nq x n_el x 2
[xs, ys] = points(geo.mesh, geo.q2.');
s = h(xs(:), ys(:));
n = numel(xs);
V = cat(3, reshape(s(1:n), size(xs)), reshape(s(n + 1:2*n), size(xs)));
end

function [lam, dir] = boundary_values(geo, k, ex, t, prm)
% hybrid vector with the Dirichlet data (face L2 projections on the
% chords) and the constrained dofs: v_hat always, u_hat unless free
n_dof = geo.n_dof;
lam = zeros(4*n_dof, 1);
dir = [];
for field = ["u", "v"]
    if field == "u" && prm.free, continue; end
    h = ex.(field);
    sp = LinearElasticityProblem(k);
    sp.Data.boundary_conditions = boundary_condition("Dirichlet", ...
        x_data = @(x, y) row(h(x(:), y(:), t), 1, size(x)), ...
        y_data = @(x, y) row(h(x(:), y(:), t), 2, size(x)));
    [x_hat, d] = get_boundary_conditions(sp, geo.mesh);
    off = 2*n_dof*(field == "v");
    lam(off + (1:2*n_dof)) = x_hat;
    dir = [dir; off + d(:)]; %#ok<AGROW>
end
end

function r = row(s, i, sz)
n = numel(s)/2;
r = reshape(s((i - 1)*n + (1:n)), sz);
end

function lam = initial_skeleton(geo, k, level, ex, t, prm)
% face L2 projection of the mean of the adjacent traces of [v; u], then
% the Dirichlet data (ex = []: zero data, a perturbation)
d2 = nchoosek(k + 2, 2);
v_hat = face_mean(geo.mesh, k, level(1:2*d2, :));
u_hat = face_mean(geo.mesh, k, level(2*d2 + 1:end, :));
lam = [reshape(u_hat(1:k + 1, :), [], 1); reshape(u_hat(k + 2:end, :), [], 1); ...
    reshape(v_hat(1:k + 1, :), [], 1); reshape(v_hat(k + 2:end, :), [], 1)];
if isempty(ex)
    is_boundary = accumarray(geo.mesh.boundaries(:), 1, [size(geo.mesh.faces, 2), 1]) == 1;
    b = reshape(repmat(is_boundary.', k + 1, 1), [], 1);
    mask = [~prm.free*b; ~prm.free*b; b; b];
    lam(logical(mask)) = 0;
else
    [data_lam, dir] = boundary_values(geo, k, ex, t, prm);
    lam(dir) = data_lam(dir);
end
end

function a_hat = face_mean(mesh, k, a)
% face L2 projection of the mean of the one-sided traces of a (2*d2 x N)
d1 = k + 1;
d2 = size(a, 1)/2;
n_el = size(mesh.elements, 2);
n_faces = size(mesh.faces, 2);
[~, ~, ~, LP, LL] = boundary_matrix(add_stabilization(mesh, 1), k);
incidence = sparse(1:3*n_el, mesh.boundaries(:), 1, 3*n_el, n_faces);
face_mass = zeros(d1, d1, 3, n_el);
for e = 1:3
    block = (e - 1)*d1 + (1:d1);
    face_mass(:, :, e, :) = reshape(LL(block, block, :), d1, d1, 1, n_el);
end
face_mass = reshape(reshape(face_mass, d1*d1, 3*n_el)*incidence, d1, d1, n_faces);
a_hat = zeros(2*d1, n_faces);
for c = 1:2
    a_c = reshape(a((c - 1)*d2 + (1:d2), :), d2, 1, n_el);
    m = reshape(pagemtimes(LP, a_c), d1, 3*n_el)*incidence;
    a_hat((c - 1)*d1 + (1:d1), :) = reshape(pagemldivide(face_mass, reshape(m, d1, 1, n_faces)), d1, n_faces);
end
end

function [level, lam, S, it, ok, history] = bdf2_step(geo, k, hp, prm, data, Xn, Xm, n, options)
% one BDF2 step: the nonlinear loop on the carve of geo (see the header);
% history = the relative change of every iteration
Z = (4*Xn - Xm)/3;
level = 2*Xn - Xm;
lam = initial_skeleton_from(geo, k, level, data);
ok = false;
history = zeros(1, 0);
for it = 1:options.picard_max_iterations
    [level_new, lam_new, S] = linear_solve(geo, k, hp, prm, data, Z, level, lam);
    change = norm([level_new(:) - level(:); lam_new - lam])/max(norm([level_new(:); lam_new]), realmin);
    history(end + 1) = change; %#ok<AGROW>
    level = level_new;
    lam = lam_new;
    if change <= options.picard_tolerance
        ok = true;
        break
    end
end
if ~ok
    warning("eulerian_hyperelastic_solid_march:picard", ...
        "Picard not converged at step %d: relative change %.2e after %d iterations.", n, change, it);
end
end

function lam = initial_skeleton_from(geo, k, level, data)
% the skeleton of an initial iterate: the face L2 projection of the mean of
% the adjacent traces of level on every face, boundary faces included,
% then the data on the constrained v_hat dofs (u_hat keeps the inner trace
% on the boundary, see the header)
d2 = nchoosek(k + 2, 2);
v_hat = face_mean(geo.mesh, k, level(1:2*d2, :));
u_hat = face_mean(geo.mesh, k, level(2*d2 + 1:end, :));
lam = [reshape(u_hat(1:k + 1, :), [], 1); reshape(u_hat(k + 2:end, :), [], 1); ...
    reshape(v_hat(1:k + 1, :), [], 1); reshape(v_hat(k + 2:end, :), [], 1)];
dir_v = data.dir(data.dir > 2*geo.n_dof);
lam(dir_v) = data.lam(dir_v);
end

function [level, lam, S] = linear_solve(geo, k, hp, prm, data, Z, level0, lam0)
% the linear solve of one iteration at the iterate (level0, lam0)
[A1, A2, Af, A3, A4, Ab] = local_operators(geo, k, hp, prm, data, Z, level0, lam0);
n_el = size(A1, 3);
Xs = pagemldivide(A1, cat(2, A2, Af));
H = pagemtimes(A3, Xs(:, 1:end - 1, :)) + A4;
r = -pagemtimes(A3, Xs(:, end, :));
nl = size(geo.idx, 1);
rows = repmat(reshape(geo.idx, nl, 1, n_el), 1, nl, 1);
cols = repmat(reshape(geo.idx, 1, nl, n_el), nl, 1, 1);
n_tot = 4*geo.n_dof;
Hg = sparse(rows(:), cols(:), H(:), n_tot, n_tot);
rg = accumarray(geo.idx(:), r(:) - Ab(:), [n_tot, 1]);
lam = data.lam;
free = true(n_tot, 1);
free(data.dir) = false;
lam(free) = Hg(free, free) \ (rg(free) - Hg(free, ~free)*lam(~free));
X = pagemtimes(Xs(:, 1:end - 1, :), reshape(lam(geo.idx), nl, 1, n_el)) + Xs(:, end, :);
X = reshape(X, [], n_el);
d2 = nchoosek(k + 2, 2);
p = geo.d2peers;
S = X(1:p, :);
level = [X(p + 3*d2 + (1:2*d2), :); X(p + (1:2*d2), :)];
end

function [A1, A2, Af, A3, A4, Ab] = local_operators(geo, k, hp, prm, data, Z, level0, lam0)
% local rows [C; K; S; M] on X = [sigma; u; gamma; v], A1 X = A2 Lam + Af,
% and the trace rows [Fu; Fv], A3 X + A4 Lam + Ab = 0, hybrid Lam =
% [u_hat; v_hat] per element
mesh = geo.mesh;
n_el = size(mesh.elements, 2);
d1 = k + 1;
d2 = nchoosek(k + 2, 2);
p = geo.d2peers;
[iS, iU, iG, iV] = deal(1:p, p + (1:2*d2), p + 2*d2 + (1:d2), p + 3*d2 + (1:2*d2));
iC = [iS, iG];
n_loc = p + 5*d2;
rho = prm.rho;
c = prm.c;
% constitutive and symmetry rows of the class
if isempty(hp)
    cls = geo.lin;
else
    n_dof = geo.n_dof;
    hp.Data.iterate = struct("displacement", level0(2*d2 + 1:end, :), "displacement_skeleton", ...
        [reshape(lam0(1:n_dof), d1, []); reshape(lam0(n_dof + 1:2*n_dof), d1, [])]);
    cls = get_matrices(hp, mesh, k);
end
% transport velocity and face weights
V0 = level0(1:2*d2, :);
Cb = convective_matrix(mesh, k, {V0(1:d2, :); V0(d2 + 1:end, :)});
CbT = permute(Cb, [2 1 3]);   % (beta.grad P_j, P_i)
lamK = reshape(lam0(geo.idx), [], 1, n_el);
bx = reshape(pagemtimes(geo.Lf, lamK(6*d1 + (1:3*d1), :, :)), [], n_el);
by = reshape(pagemtimes(geo.Lf, lamK(9*d1 + (1:3*d1), :, :)), [], n_el);
bn = geo.nx.*bx + geo.ny.*by;
a_v = abs(bn);
a_u = (bn + sqrt(bn.^2 + prm.delta^2))/2 + prm.tau_u;
[PPv, LPv] = weighted_face(geo, prm.tau_v + rho*(a_v - bn));
[~, LPv2, LLv2] = weighted_face(geo, prm.tau_v + rho*a_v);
[PPu, LPu] = weighted_face(geo, a_u - bn);
w_Fu = a_u;
Ab = zeros(12*d1, 1, n_el);
kinematic = prm.char && isfield(data, "kinematic");
if kinematic
    % g_u from the iterate: the BDF2 step of the kinematic equation at the
    % face points, Z_u + (2 dt/3)(v - (grad u) v) (data.kinematic.Zu, .dt)
    ck = 2*data.kinematic.dt/3;
    U0 = level0(2*d2 + 1:end, :);
    Vf = {geo.Pf*V0(1:d2, :), geo.Pf*V0(d2 + 1:end, :)};
    Gf = cell(2, 2);   % {c, d}: d_d u_c at the face points
    for ci = 1:2
        Uc = reshape(U0((ci - 1)*d2 + (1:d2), :), d2, 1, n_el);
        Gf{ci, 1} = reshape(pagemtimes(geo.Pfx, Uc), [], n_el);
        Gf{ci, 2} = reshape(pagemtimes(geo.Pfy, Uc), [], n_el);
        data.gu_face(:, :, ci) = geo.Pf*data.kinematic.Zu((ci - 1)*d2 + (1:d2), :) ...
            + ck*(Vf{ci} - Gf{ci, 1}.*Vf{1} - Gf{ci, 2}.*Vf{2});
    end
end
if prm.char
    % <a_u (u - u_hat) + a_in (g_u - u_hat), mu> on boundary faces
    a_in = (-bn + sqrt(bn.^2 + prm.delta^2))/2.*(geo.bnd == 1);
    [~, LPu2] = weighted_face(geo, w_Fu);
    [~, ~, LLu2] = weighted_face(geo, w_Fu + a_in);
    W = geo.wa.*a_in;
    for comp = 1:2
        Ab((comp - 1)*3*d1 + (1:3*d1), 1, :) = pagemtimes(geo.Lf, "transpose", ...
            reshape(W.*data.gu_face(:, :, comp), [], 1, n_el), "none");
    end
else
    if prm.free, w_Fu(geo.bnd == 1) = 1; end
    [~, LPu2, LLu2] = weighted_face(geo, w_Fu);
end
M = geo.M;

A1 = zeros(n_loc, n_loc, n_el);
A2 = zeros(n_loc, 12*d1, n_el);
Af = zeros(n_loc, 1, n_el);
A1(iC, [iS, iU, iG], :) = cls{1}(iC, :, :);
A2(iC, 1:6*d1, :) = cls{2}(iC, :, :);
Af(iC, :, :) = cls{3}(iC, :, :);
% (K)
A1(iU, iU, :) = blk2(c*M + CbT + PPu);
A1(iU, iV, :) = -blk2(M);
A2(iU, 1:6*d1, :) = blk2(permute(LPu, [2 1 3]));
Af(iU, :, :) = c*pagemtimes(blk2(M), reshape(Z(2*d2 + 1:end, :), [], 1, n_el)) + data.g;
% (M)
A1(iV, iS, :) = -permute(geo.tU, [2 1 3]);
A1(iV, iV, :) = blk2(c*rho*M + rho*CbT + PPv);
A2(iV, 6*d1 + 1:end, :) = blk2(permute(LPv, [2 1 3]));
Af(iV, :, :) = c*rho*pagemtimes(blk2(M), reshape(Z(1:2*d2, :), [], 1, n_el)) + data.f;
% trace rows: (Fu) on u_hat, (Fv) on v_hat
A3 = zeros(12*d1, n_loc, n_el);
A4 = zeros(12*d1, 12*d1, n_el);
A3(1:6*d1, iU, :) = blk2(LPu2);
A4(1:6*d1, 1:6*d1, :) = -blk2(LLu2);
A3(6*d1 + 1:end, iS, :) = geo.Sn;
A3(6*d1 + 1:end, iV, :) = -blk2(LPv2);
A4(6*d1 + 1:end, 6*d1 + 1:end, :) = blk2(LLv2);
if ~prm.newton, return, end
% Newton: the derivatives of the transport terms with respect to beta = v
% and beta_hat = v_hat (through beta_hat.n in the weights), at the iterate:
% JX on the element v, JL on the hybrid v_hat in the local rows, TL in the
% trace rows; the system for the new iterate is (A + J - A) x = b + (J - A) x0
U0 = level0(2*d2 + 1:end, :);
trace_hat = @(off) reshape(pagemtimes(geo.Lf, lamK(off + (1:3*d1), :, :)), [], n_el);
du = {geo.Pf*U0(1:d2, :) - trace_hat(0), geo.Pf*U0(d2 + 1:end, :) - trace_hat(3*d1)};   % u0 - u_hat0
dv = {geo.Pf*V0(1:d2, :) - bx, geo.Pf*V0(d2 + 1:end, :) - by};                          % v0 - v_hat0
nn = {geo.nx, geo.ny};
r = sqrt(bn.^2 + prm.delta^2);
da_u = (1 + bn./r)/2;
s_v = sign(bn);
w_Fu = da_u;
if prm.free, w_Fu(geo.bnd == 1) = 0; end
Gu = gradient_mass(geo, k, U0);
Gv = gradient_mass(geo, k, V0);
JX = zeros(n_loc, n_loc, n_el);
JL = zeros(n_loc, 12*d1, n_el);
TL = zeros(12*d1, 12*d1, n_el);
for ci = 1:2
    rU = iU((ci - 1)*d2 + (1:d2));
    rV = iV((ci - 1)*d2 + (1:d2));
    rFu = (ci - 1)*3*d1 + (1:3*d1);
    rFv = 6*d1 + rFu;
    w_char = 0;
    if prm.char
        w_char = (-1 + bn./r)/2.*(data.gu_face(:, :, ci) - trace_hat((ci - 1)*3*d1)).*(geo.bnd == 1);
    end
    for di = 1:2
        cV = iV((di - 1)*d2 + (1:d2));
        cL = 6*d1 + (di - 1)*3*d1 + (1:3*d1);
        JX(rU, cV, :) = Gu{ci, di};
        JX(rV, cV, :) = rho*Gv{ci, di};
        [~, LP] = weighted_face(geo, (da_u - 1).*nn{di}.*du{ci});
        JL(rU, cL, :) = permute(LP, [2 1 3]);
        [~, LP] = weighted_face(geo, rho*(s_v - 1).*nn{di}.*dv{ci});
        JL(rV, cL, :) = permute(LP, [2 1 3]);
        [~, ~, LL] = weighted_face(geo, (w_Fu.*du{ci} + w_char).*nn{di});
        TL(rFu, cL, :) = LL;
        [~, ~, LL] = weighted_face(geo, -rho*s_v.*nn{di}.*dv{ci});
        TL(rFv, cL, :) = LL;
    end
end
A1 = A1 + JX;
A2 = A2 - JL;
Af = Af + pagemtimes(JX(:, iV, :), reshape(V0, [], 1, n_el)) + pagemtimes(JL, lamK);
A4 = A4 + TL;
Ab = Ab - pagemtimes(TL, lamK);
if kinematic
    % the derivative of <a_in g_u, mu> with respect to the element v and u:
    % d g_u,c = ck (dv_c - (grad u_c).dv - (grad du_c).v)
    TX = zeros(12*d1, n_loc, n_el);
    for ci = 1:2
        rFu = (ci - 1)*3*d1 + (1:3*d1);
        for di = 1:2
            [~, LP] = weighted_face(geo, a_in.*ck.*((ci == di) - Gf{ci, di}));
            TX(rFu, iV((di - 1)*d2 + (1:d2)), :) = LP;
        end
        W = reshape(-ck*geo.wa.*a_in, [], 1, n_el);
        TX(rFu, iU((ci - 1)*d2 + (1:d2)), :) = pagemtimes(geo.Lf, "transpose", ...
            W.*(reshape(Vf{1}, [], 1, n_el).*geo.Pfx + reshape(Vf{2}, [], 1, n_el).*geo.Pfy), "none");
    end
    A3 = A3 + TX;
    Ab = Ab - pagemtimes(TX, reshape([zeros(p, n_el); U0; zeros(d2, n_el); V0], [], 1, n_el));
end
end

function G = gradient_mass(geo, k, U)
% {c, d}: (d_d U_c P_j, P_i)_K of the element field U (2*d2 x n_el),
% d2 x d2 x n_el each
mesh = geo.mesh;
n_el = size(mesh.elements, 2);
d2 = nchoosek(k + 2, 2);
[rs, w] = rule(ceil(3*k/2) + 1);
[P, Pr, Ps] = dubiner_basis(k, rs(1, :), rs(2, :));
F = get_affine_map(mesh, output = "array");
dt = F(1, 1, :).*F(2, 2, :) - F(1, 2, :).*F(2, 1, :);
Pd = {(F(2, 2, :).*Pr - F(2, 1, :).*Ps)./dt, (-F(1, 2, :).*Pr + F(1, 1, :).*Ps)./dt};
WP = (w.*reshape(mesh.volume, 1, 1, [])).*P;
G = cell(2, 2);
for c = 1:2
    Uc = reshape(U((c - 1)*d2 + (1:d2), :), d2, 1, n_el);
    for d = 1:2
        g = pagemtimes(Pd{d}, Uc);
        G{c, d} = pagemtimes(P, "transpose", g.*WP, "none");
    end
end
end

function err = jacobian_check(geo, k, hp, prm, data, Z, level0, lam0)
% relative difference between the central difference of the residual
% (A(x) x - b(x), the free rows) in a random direction d and the Newton
% operator applied to d, for the steps 1e-3 ... 1e-8 relative to |x0|
n_el = size(geo.mesh.elements, 2);
d2 = nchoosek(k + 2, 2);
p = geo.d2peers;
n_loc = p + 5*d2;
[iU, iV] = deal(p + (1:2*d2), p + 3*d2 + (1:2*d2));
stream = RandStream("mt19937ar", Seed = 1);
X0 = 0.1*randn(stream, n_loc, n_el);
X0(iU, :) = level0(2*d2 + 1:end, :);
X0(iV, :) = level0(1:2*d2, :);
dX = randn(stream, n_loc, n_el);
dL = randn(stream, size(lam0));
dL(data.dir) = 0;
s = norm([X0(:); lam0])/norm([dX(:); dL]);
[dX, dL] = deal(s*dX, s*dL);
pn = prm;
pn.newton = true;
[A1, A2, ~, A3, A4] = local_operators(geo, k, hp, pn, data, Z, level0, lam0);
Jd = residual_rows(geo, data, A1, A2, 0, A3, A4, 0, dX, dL);
pr = prm;
pr.newton = false;
steps = 10.^(-3:-1:-8);
err = zeros(size(steps));
for i = 1:numel(steps)
    e = steps(i);
    fd = (residual(X0 + e*dX, lam0 + e*dL) - residual(X0 - e*dX, lam0 - e*dL))/(2*e);
    err(i) = norm(fd - Jd)/norm(Jd);
end
    function r = residual(X, L)
        [B1, B2, Bf, B3, B4, Bb] = local_operators(geo, k, hp, pr, data, Z, [X(iV, :); X(iU, :)], L);
        r = residual_rows(geo, data, B1, B2, Bf, B3, B4, Bb, X, L);
    end
end

function R = residual_rows(geo, data, A1, A2, Af, A3, A4, Ab, X, lam)
% [A1 X - A2 Lam - Af (element rows); sum_K (A3 X + A4 Lam + Ab) (free
% trace rows)]
n_el = size(geo.mesh.elements, 2);
X = reshape(X, [], 1, n_el);
LK = reshape(lam(geo.idx), [], 1, n_el);
rl = pagemtimes(A1, X) - pagemtimes(A2, LK) - Af;
rt = pagemtimes(A3, X) + pagemtimes(A4, LK) + Ab;
n_tot = 4*geo.n_dof;
rg = accumarray(geo.idx(:), rt(:), [n_tot, 1]);
free = true(n_tot, 1);
free(data.dir) = false;
R = [rl(:); rg(free)];
end

function Pp = peers_points(mesh, k, rs)
% PEERS basis at the reference points rs of every element, as
% (nq*4) x nP x n_el, rows q fastest then the components 11, 21, 12, 22
G = get_affine_map(mesh, inverse = true);
P = PEERS_basis(k, rs(1, :), rs(2, :), G);
[nq, nP] = size(P, 1:2);
Pp = reshape(permute(P, [1 3 4 2 5]), nq*4, nP, []);
end

function vals = peers_values(mesh, k, S, rs)
Pp = peers_points(mesh, k, rs);
vals = reshape(pagemtimes(Pp, reshape(S, size(S, 1), 1, [])), size(Pp, 1), []);
end

function L = lifting(geo, k, U, lam)
% PEERS dofs of L_h(u, u_hat)
n_el = size(U, 2);
uhK = lam(geo.idx(1:6*(k + 1), :));
L = pagemtimes(geo.R_u, reshape(U, [], 1, n_el)) + pagemtimes(geo.R_hat, reshape(uhK, [], 1, n_el));
L = reshape(L, [], n_el);
end

function e = strain_energy(H, prm)
% 1/2 (C eps(H), eps(H)) pointwise, H values (nq*4) x n (11, 21, 12, 22)
nq = size(H, 1)/4;
h = @(i) H((i - 1)*nq + (1:nq), :);
s12 = (h(2) + h(3))/2;
tr = h(1) + h(4);
e = prm.mu*(h(1).^2 + h(4).^2 + 2*s12.^2) + prm.lambda/2*tr.^2;
end

function E = linearized_energy(geo, k, prm, level, lam)
% 1/2 rho |v|^2 + 1/2 (C eps(L_h), eps(L_h)) of (v, u, u_hat)
E = deviation(geo, k, prm, level, lam, [], 0);
end

function E = deviation(geo, k, prm, level, lam, ex, t)
% 1/2 rho |v_h - v|^2 + 1/2 (C eps(L_h - grad u), eps(L_h - grad u))
% (ex = []: the energy of (v_h, L_h) itself)
mesh = geo.mesh;
d2 = nchoosek(k + 2, 2);
[rs, w] = rule(k + 3);
nq = numel(w);
D = dubiner_basis(k, rs);
vh = [D*level(1:d2, :); D*level(d2 + (1:d2), :)];
Lh = peers_values(mesh, k, lifting(geo, k, level(2*d2 + 1:end, :), lam), rs);
if ~isempty(ex)
    [xs, ys] = points(mesh, rs);
    vh = vh - unstack(ex.v(xs(:), ys(:), t), xs);
    Lh = Lh - unstack(ex.grad_u(xs(:), ys(:), t), xs);
end
W = w.*mesh.volume(:).';
ev = sum(W.*(vh(1:nq, :).^2 + vh(nq + 1:end, :).^2), "all");
es = sum(W.*strain_energy(Lh, prm), "all");
E = prm.rho*ev/2 + es;
end

function V = unstack(s, xs)
% values stacked [component 1 (numel(xs)); component 2; ...] at the points
% xs (nq x n) -> (nq*n_c) x n, rows q fastest then component
[nq, n] = size(xs);
V = reshape(permute(reshape(s, nq, n, []), [1 3 2]), [], n);
end

function [e, nrm] = errors(geo, k, level, S, ex, t)
% L2 errors and exact norms [v, u, sigma] at time t
mesh = geo.mesh;
d2 = nchoosek(k + 2, 2);
[rs, w] = rule(k + 3);
[xs, ys] = points(mesh, rs);
D = dubiner_basis(k, rs);
W = w.*mesh.volume(:).';
nq = numel(w);
sq = @(a, nc) reshape(sum(reshape(a, nq, nc, []).^2, 2), nq, []);
num = {kron(eye(2), D)*level(1:2*d2, :), kron(eye(2), D)*level(2*d2 + 1:end, :), ...
    peers_values(mesh, k, S, rs)};
exact = {unstack(ex.v(xs(:), ys(:), t), xs), unstack(ex.u(xs(:), ys(:), t), xs), ...
    unstack(ex.sigma(xs(:), ys(:), t), xs)};
nc = [2, 2, 4];
[e, nrm] = deal(zeros(1, 3));
for i = 1:3
    e(i) = sqrt(sum(W.*sq(num{i} - exact{i}, nc(i)), "all"));
    nrm(i) = sqrt(sum(W.*sq(exact{i}, nc(i)), "all"));
end
end

function b = field_bounds(mesh, k, ex, t, c_s)
% [max |grad u| (Frobenius), max |v|/c_s, min det(I - grad u)] at the
% quadrature points of mesh
[rs, ~] = rule(k + 3);
[xs, ys] = points(mesh, rs);
n = numel(xs);
G = reshape(ex.grad_u(xs(:), ys(:), t), n, 4);
v = reshape(ex.v(xs(:), ys(:), t), n, 2);
b = [max(sqrt(sum(G.^2, 2))), max(sqrt(sum(v.^2, 2)))/c_s, ...
    min((1 - G(:, 1)).*(1 - G(:, 4)) - G(:, 2).*G(:, 3))];
end
