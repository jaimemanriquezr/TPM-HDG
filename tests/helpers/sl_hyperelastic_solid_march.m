function out = sl_hyperelastic_solid_march(k, N, amplitude, dt, n_steps, options)
%SL_HYPERELASTIC_SOLID_MARCH  Semi-Lagrangian Eulerian hyperelastic solid on a moving snapped carve.
%
%   out = sl_hyperelastic_solid_march(k, N, A, dt, n_steps, rotation = 0, ...
%       law = "neo_hookean", fields = "standard", scheme = "bdf2", ...)
%
%   The solid alone, fully Eulerian, in the current configuration Omega(t):
%   the disc of the given radius and centre c(t) = c0 + A sin(2 pi t) e_x,
%   carved at every t_n from the background StructuredSimplicialGrid(N, N,
%   "crossed") on [0, 1]^2 with its interface vertices snapped onto the
%   circle (the carve of eulerian_hyperelastic_solid_march); every boundary
%   face is Dirichlet. With D_t = d_t + v.grad and G = grad u (x - u = X
%   the reference map):
%     rho D_t v - div sigma = f,   D_t u = v + g,   sigma = sigma_E(G),
%   so that D_t G = (I - G) grad v + grad g; law = "neo_hookean": the
%   Cauchy stress of HyperelasticProblem.first_piola (material
%   "neo_hookean_eulerian"); law = "linear": sigma = C eps(G). Split law
%   sigma_E = C eps + N_E, N_E symmetric.
%
%   State (degree k): v in P_k^2 and L in S_h (the PEERS space, full 2 x 2)
%   per element, plus the passive displacement u_P in P_k^2, which is
%   carried and reported but does not enter the dynamics. The stress never
%   comes from a lifting of a displacement: L is transported and updated by
%   the PEERS lifting of (v, v_hat), v_hat the trace of the momentum flux
%   balance.
%
%   Transport (self-transport). Feet x in Omega_n+1 -> y = phi(x) in
%   Omega_n by the trapezoid rule y = x - dt/2 (beta(x) + v^n(y)), beta =
%   2 v^n - v^n-1 (both carried to Omega_n+1 by the zero-shift
%   LinearElasticityProblem.translated_field_between), grad phi =
%   (I + dt/2 grad v^n(y))^-1 (I - dt/2 grad beta(x)). Transported fields
%   v~ = Pi[v o phi], u_P~ = Pi[u_P o phi], L~ = Pi_S[(L o phi) grad phi]
%   (the gradient of the transported displacement). On every new element
%   K, the affine interpolant phi_K of the feet at K's vertices decides
%   which old element owns which part of K: phi_K(K) is clipped against the
%   old triangles, the pieces are mapped back to K and integrated with a
%   Gauss rule whose integrand uses the true foot phi(x) and the polynomial
%   of the owning old element (also for v^n(y) in the foot iteration).
%   Where phi_K(K) leaves the old carve, the polynomial of one old element
%   E is extrapolated (E = largest overlap, else the nearest centroid). L
%   travels in its exact P_k+1 Dubiner form and is projected onto S_h on K.
%   The second BDF2 level is the previous step's transported level,
%   transported once more.
%
%   Step on the carve of t_n+1 (scheme = "bdf2": c = 3/(2 dt), H = (4 X~^n
%   - X~^n-1)/3; scheme = "cn": c = 2/dt, H = X~ + R~/c, R~ the transported
%   residual R^n = c (X^n - X~) - R~ of the previous step, L-type with the
%   weight grad phi), tau in S_h, eta skew P_k, w in P_k^2:
%    (C) (A sigma, tau) + (gamma, tau) - (L, tau) - (A N_E(L), tau) = 0
%    (S) (sigma, eta) = 0
%    (L) c (L, tau) + (v, div tau) - <v_hat, tau n> = c (H_L, tau) + (grad g, tau)
%    (M) c rho (v, w) - (div sigma, w) + <tau_v (v - v_hat), w> = c rho (H_v, w) + (f, w)
%    (F) sum_K <sigma n - tau_v (v - v_hat), mu> = 0 on interior faces,
%        v_hat = P g_v (face L2 projection of the exact velocity) on the chords;
%   u_P = (4 u~^n - u~^n-1)/3 + 2 dt/3 (v + g) (bdf2) or u_P = u~ + dt/2
%   (v + g + (v + g)~) (cn); f, g at t_n+1; tau_v = rho c_s by default.
%   Only (C) is nonlinear: Newton on N_E at L (first_piola's derivative),
%   one statically condensed solve per iteration, stopped at a relative
%   change of (sigma, L, v, v_hat) <= newton_tolerance.
%
%   Exact fields: "standard" and "quadratic" as in
%   eulerian_hyperelastic_solid_march; "patch": a uniform material velocity
%   v = patch_velocity and u = U(x - v t, t), U of degree patch_degree in x
%   and (1 + t^2) in t, which the scheme reproduces (feet, transport, BDF2
%   and CN are exact on it; the rotation does not enter). Start: the L2
%   projections of the exact (v, grad u, u) at t = 0 on carve(0) and at
%   t = -dt on carve(-dt) (scheme = "cn": also the exact rates (D_t v,
%   grad v + grad g) at t = 0 as R^0).
%   transport_check = omega ~= 0: no march; the patch level at t = 0 is
%   transported once to carve(dt) with the affine velocity v = patch_velocity
%   + omega (x - c0)^perp (trapezoid feet = an affine Cayley map), and
%   out.transport_error is the largest relative L2 difference of (v, L, u)
%   from the projections of the exactly composed fields.
%   initial_perturbation = [seed, amplitude]: a twin state with its own
%   feet, its level 0 (and the transported level -1) plus amplitude times a
%   random (dv, dL) of unit energy 1/2 rho |dv|^2 + 1/2 (C eps(dL), eps(dL)).
%
%   out: e_v, e_L, e_sigma, e_u (L2 errors at T: v, L against grad u, the
%   S_h stress against the exact Cauchy stress, u_P), n_v, n_L, n_sigma,
%   n_u (exact norms); step_errors (record_errors); iterations, converged
%   (n_steps x 1, both states); uncovered (n_steps x 2); energy (n_steps +
%   1, 1/2 rho |v|^2 + 1/2 (C eps(L), eps(L))); perturbation_energy (the
%   same of the twin difference); field_bounds = [max |grad u|, max |v|/c_s,
%   min det(I - grad u)] of the exact fields; max_L (largest |L| at the
%   quadrature points over the run); final v, L, sigma, u, mesh; seconds.
%
%See also: eulerian_hyperelastic_solid_march, velocity_stress_solid_march,
%   HyperelasticProblem.first_piola
arguments
    k (1,1) double
    N (1,1) double
    amplitude (1,1) double
    dt (1,1) double
    n_steps (1,1) double
    options.rotation (1,1) double = 0
    options.law (1,1) string {mustBeMember(options.law, ["neo_hookean", "linear"])} = "neo_hookean"
    options.fields (1,1) string {mustBeMember(options.fields, ["standard", "quadratic", "patch"])} = "standard"
    options.scheme (1,1) string {mustBeMember(options.scheme, ["bdf2", "cn"])} = "bdf2"
    options.patch_degree double = []
    options.patch_velocity (2,1) double = [0.3; -0.2]
    options.w (2,1) string = ["0.02*sin(5*X + 3*Y)*cos(3*t)"; "0.02*cos(4*X - 5*Y)*(1 + t^2)"]
    options.tau_v double = []
    options.newton_tolerance (1,1) double = 1e-11
    options.newton_max_iterations (1,1) double = 20
    options.initial_perturbation double = []
    options.transport_check (1,1) double = 0
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
[lambda, mu] = LinearElasticityProblem.get_lame_parameters(options.E, options.nu);
prm = struct("rho", rho, "lambda", lambda, "mu", mu, "law", options.law, ...
    "tau_v", default_to(options.tau_v, rho*c_s), "dt", dt, ...
    "c", (options.scheme == "bdf2")*3/(2*dt) + (options.scheme == "cn")*2/dt, ...
    "tol", options.newton_tolerance, "maxit", options.newton_max_iterations);
if isempty(options.patch_degree), options.patch_degree = k; end
cen = @(t) options.center + [amplitude*sin(2*pi*t); 0];
moving = amplitude ~= 0;
carve = @(t) carve_solid(local_background(N, cen(t), options.radius), ...
    CircularCurve(cen(t), options.radius), options.snap);
ex = exact_fields(options, amplitude, lambda, mu);
sp0 = LinearElasticityProblem(k);
to_new = @(old, new, U) sp0.translated_field_between(old, new, k, k, U, [0; 0]);

geo = geometry(carve(0), k, options);
if options.transport_check ~= 0
    out = transport_check(geo, geometry(carve(dt), k, options), k, dt, ex, options);
    out.seconds = toc(wall);
    return
end
geo_prev = geo;
if moving, geo_prev = geometry(carve(-dt), k, options); end
st = initial_level(geo, k, ex, 0);
stm = initial_level(geo_prev, k, ex, -dt);
Y = transport_state(transport_setup(geo_prev, geo, k, dt, st.v, stm.v), stm, k);
v_prev = stm.v;
R = [];
if options.scheme == "cn"
    R = struct("v", project_vector(geo.mesh, k, @(x, y) ex.dtv(x, y, 0)), ...
        "L", project_peers(geo, k, @(x, y) ex.grad_vg(x, y, 0)));
end
twin = ~isempty(options.initial_perturbation);
if twin
    stream = RandStream("mt19937ar", Seed = options.initial_perturbation(1));
    d = struct("v", randn(stream, size(st.v)), "L", randn(stream, size(st.L)));
    d0 = options.initial_perturbation(2)/sqrt(energy(geo, k, prm, d.v, d.L));
    st2 = st;
    st2.v = st.v + d0*d.v;
    st2.L = st.L + d0*d.L;
    Y2 = Y;
    Y2.v = Y.v + d0*d.v;
    Y2.L = Y.L + d0*d.L;
    v_prev2 = v_prev;
    if moving
        v_prev2 = v_prev + d0*to_new(geo.mesh, geo_prev.mesh, d.v);
    else
        v_prev2 = v_prev + d0*d.v;
    end
    R2 = R;
end

energies = nan(n_steps + 1, 1);
perturbation_energy = nan(n_steps + 1, 1);
energies(1) = energy(geo, k, prm, st.v, st.L);
if twin, perturbation_energy(1) = energy(geo, k, prm, st2.v - st.v, st2.L - st.L); end
bounds = field_bounds(geo.mesh, k, ex, 0, c_s);
iterations = zeros(n_steps, 1);
converged = true(n_steps, 1);
uncovered = zeros(n_steps, 2);
step_errors = nan(n_steps, 4);
max_L = 0;
for n = 1:n_steps
    t1 = n*dt;
    geo_new = geo;
    if moving, geo_new = geometry(carve(t1), k, options); end
    data = step_data(geo_new, k, ex, t1);
    [st_new, Y_new, R_new, S, it, ok, info] = sl_step(geo_prev, geo, geo_new, k, prm, data, st, Y, ...
        v_prev, R, moving, to_new, options.scheme);
    iterations(n) = it;
    converged(n) = ok;
    uncovered(n, :) = [max(info.uncovered), sum(info.uncovered(:).*geo_new.mesh.volume(:))];
    if twin
        [st2_new, Y2, R2, ~, ~, ok2] = sl_step(geo_prev, geo, geo_new, k, prm, data, st2, Y2, ...
            v_prev2, R2, moving, to_new, options.scheme);
        converged(n) = ok && ok2;
        v_prev2 = st2.v;
        st2 = st2_new;
        perturbation_energy(n + 1) = energy(geo_new, k, prm, st2.v - st_new.v, st2.L - st_new.L);
    end
    v_prev = st.v;
    [st, Y, R] = deal(st_new, Y_new, R_new);
    [geo_prev, geo] = deal(geo, geo_new);
    energies(n + 1) = energy(geo, k, prm, st.v, st.L);
    b = field_bounds(geo.mesh, k, ex, t1, c_s);
    bounds = [max(bounds(1:2), b(1:2)), min(bounds(3), b(3))];
    Lq = pagemtimes(geo.Pp, reshape(st.L, [], 1, size(st.L, 2)));
    nq = size(geo.Pp, 1)/4;
    max_L = max(max_L, max(sqrt(sum(reshape(Lq, nq, 4, []).^2, 2)), [], "all"));
    if options.record_errors
        [e, ~] = errors(geo, k, st, S, ex, t1);
        step_errors(n, :) = e;
    end
    if ~ok, break, end
end
T = n*dt;
[e, nrm] = errors(geo, k, st, S, ex, T);
out = struct("e_v", e(1), "e_L", e(2), "e_sigma", e(3), "e_u", e(4), "n_v", nrm(1), ...
    "n_L", nrm(2), "n_sigma", nrm(3), "n_u", nrm(4), "step_errors", step_errors, ...
    "iterations", iterations, "converged", converged, "uncovered", uncovered, ...
    "energy", energies, "perturbation_energy", perturbation_energy, "field_bounds", bounds, ...
    "max_L", max_L, "v", st.v, "L", st.L, "sigma", S, "u", st.u, "mesh", geo.mesh, ...
    "n_elements", size(geo.mesh.elements, 2), "steps_done", n, "seconds", toc(wall));
end

function x = default_to(x, value)
if isempty(x), x = value; end
end

%% ------------------------------------------------------------------ step
function [st_new, Y_new, R_new, S, it, ok, info] = sl_step(geo_prev, geo, geo_new, k, prm, data, ...
    st, Y, v_prev, R, moving, to_new, scheme)
% one step from the carve of geo (state st, transported level Y of the
% previous step, v_prev on geo_prev) to the carve of geo_new
if moving
    beta = 2*to_new(geo.mesh, geo_new.mesh, st.v) - to_new(geo_prev.mesh, geo_new.mesh, v_prev);
else
    beta = 2*st.v - v_prev;
end
tr = transport_setup(geo, geo_new, k, prm.dt, beta, st.v);
info.uncovered = tr.uncovered;
Xn = transport_state(tr, st, k);
c = prm.c;
if scheme == "bdf2"
    Xm = transport_state(tr, Y, k);
    H = struct("v", (4*Xn.v - Xm.v)/3, "L", (4*Xn.L - Xm.L)/3);
    L0 = 2*Xn.L - Xm.L;
else
    Rt = transport_state(tr, R, k);
    H = struct("v", Xn.v + Rt.v/c, "L", Xn.L + Rt.L/c);
    L0 = H.L;
end
[v, L, S, it, ok] = step_solve(geo_new, k, prm, data, H, L0);
st_new = struct("v", v, "L", L);
R_new = [];
if scheme == "bdf2"
    st_new.u = (4*Xn.u - Xm.u)/3 + 2*prm.dt/3*(v + data.g);
else
    st_new.u = Xn.u + prm.dt/2*(v + data.g + Xn.v + Xn.g);
    R_new = struct("v", c*(v - Xn.v) - Rt.v, "L", c*(L - Xn.L) - Rt.L);
end
st_new.g = data.g;
Y_new = Xn;
end

function [v, L, S, it, ok] = step_solve(geo, k, prm, data, H, L0)
% Newton on (C); the other rows are linear (see the header)
mesh = geo.mesh;
n_el = size(mesh.elements, 2);
d1 = k + 1;
d2 = nchoosek(k + 2, 2);
p = geo.p;
[iS, iG, iL, iV] = deal(1:p, p + (1:d2), p + d2 + (1:p), 2*p + d2 + (1:2*d2));
n_loc = 2*p + 3*d2;
c = prm.c;
rho = prm.rho;
tv = prm.tau_v;
A1 = zeros(n_loc, n_loc, n_el);
A2 = zeros(n_loc, 6*d1, n_el);
Af = zeros(n_loc, 1, n_el);
A1(iS, iS, :) = geo.tS;
A1(iS, iG, :) = geo.tG;
A1(iG, iS, :) = permute(geo.tG, [2 1 3]);
A1(iL, iL, :) = c*geo.MP;
A1(iL, iV, :) = geo.tU;
A2(iL, :, :) = geo.tH;
Af(iL, :, :) = c*pagemtimes(geo.MP, reshape(H.L, [], 1, n_el)) + data.gradg;
A1(iV, iS, :) = -permute(geo.tU, [2 1 3]);
A1(iV, iV, :) = blk2(c*rho*geo.M + tv*geo.tPP);
A2(iV, :, :) = tv*blk2(permute(geo.tLP, [2 1 3]));
Af(iV, :, :) = c*rho*pagemtimes(blk2(geo.M), reshape(H.v, [], 1, n_el)) + data.f;
A3 = zeros(6*d1, n_loc, n_el);
A3(:, iS, :) = geo.Sn;
A3(:, iV, :) = -tv*blk2(geo.tLP);
A4 = tv*blk2(geo.tLL);
lam = data.lam;
L = L0;
X = [];
ok = false;
for it = 1:prm.maxit
    A1i = A1;
    Afi = Af;
    if prm.law == "linear"
        A1i(iS, iL, :) = -geo.MP;
    else
        [J, r] = remainder_terms(geo, prm, L);
        A1i(iS, iL, :) = -geo.MP - J;
        Afi(iS, :, :) = r - pagemtimes(J, reshape(L, [], 1, n_el));
    end
    Xs = pagemldivide(A1i, cat(2, A2, Afi));
    Hm = pagemtimes(A3, Xs(:, 1:end - 1, :)) + A4;
    rr = -pagemtimes(A3, Xs(:, end, :));
    nl = size(geo.idx, 1);
    rows = repmat(reshape(geo.idx, nl, 1, n_el), 1, nl, 1);
    cols = repmat(reshape(geo.idx, 1, nl, n_el), nl, 1, 1);
    n_tot = 2*geo.n_dof;
    Hg = sparse(rows(:), cols(:), Hm(:), n_tot, n_tot);
    rg = accumarray(geo.idx(:), rr(:), [n_tot, 1]);
    lam_new = data.lam;
    free = true(n_tot, 1);
    free(data.dir) = false;
    lam_new(free) = Hg(free, free) \ (rg(free) - Hg(free, ~free)*lam_new(~free));
    X_new = reshape(pagemtimes(Xs(:, 1:end - 1, :), reshape(lam_new(geo.idx), nl, 1, n_el)) ...
        + Xs(:, end, :), n_loc, n_el);
    if prm.law == "linear"
        X = X_new;
        ok = true;
        break
    end
    if ~isempty(X)
        change = norm([X_new(:) - X(:); lam_new - lam])/max(norm([X_new(:); lam_new]), realmin);
    else
        change = inf;
    end
    X = X_new;
    lam = lam_new;
    L = X(iL, :);
    if change <= prm.tol
        ok = true;
        break
    end
end
if ~ok
    warning("sl_hyperelastic_solid_march:newton", "Newton not converged after %d iterations.", it);
end
S = X(iS, :);
L = X(iL, :);
v = X(iV, :);
end

function [J, r] = remainder_terms(geo, prm, L)
% r = (A N_E(L), tau) and its derivative J = (A DN_E(L)[tau_j], tau_i)
n_el = size(L, 2);
nq = size(geo.Pp, 1)/4;
Lq = reshape(pagemtimes(geo.Pp, reshape(L, [], 1, n_el)), nq, 4, n_el);
G = reshape(permute(Lq, [1 3 2]), [], 4);
[s, ds] = HyperelasticProblem.first_piola(G, prm.lambda, prm.mu, "neo_hookean_eulerian");
[Cm, Am] = elastic_matrices(prm);
AN = (s - G*Cm.')*Am.';                                          % A N_E, n x 4
ADN = zeros(size(G, 1), 4, 4);
for a = 1:4
    for b = 1:4
        ADN(:, a, b) = sum(Am(a, :).*(reshape(ds(:, :, b), [], 4) - Cm(:, b).'), 2);
    end
end
W = reshape(geo.wv, nq, 1, n_el);
P4 = reshape(geo.Pp, nq, 4, [], n_el);                           % q, comp, basis, el
AN = permute(reshape(AN, nq, n_el, 4), [1 3 2]);                 % q, comp, el
r = zeros(geo.p, 1, n_el);
for a = 1:4
    r = r + pagemtimes(reshape(P4(:, a, :, :), nq, [], n_el), "transpose", W.*AN(:, a, :), "none");
end
ADN = permute(reshape(ADN, nq, n_el, 4, 4), [1 3 4 2]);          % q, a, b, el
J = zeros(geo.p, geo.p, n_el);
for b = 1:4
    for a = 1:4
        Pa = reshape(P4(:, a, :, :), nq, [], n_el);
        Pb = reshape(P4(:, b, :, :), nq, [], n_el);
        J = J + pagemtimes(Pa, "transpose", W.*reshape(ADN(:, a, b, :), nq, 1, n_el).*Pb, "none");
    end
end
end

function [Cm, Am] = elastic_matrices(prm)
% C eps(H) = Cm H and A s = Am s on the column-major 4-vectors [11 21 12 22]
I4 = eye(4);
Cm = zeros(4);
Am = zeros(4);
c1 = 1/(2*prm.mu);
c2 = prm.lambda/(4*prm.mu*(prm.lambda + prm.mu));
for b = 1:4
    H = I4(:, b);
    tr = H(1) + H(4);
    Cm(:, b) = prm.mu*(H + H([1 3 2 4])) + prm.lambda*tr*[1; 0; 0; 1];
    Am(:, b) = c1*H - c2*tr*[1; 0; 0; 1];
end
end

%% ------------------------------------------------------------- transport
function tr = transport_setup(geo_old, geo_new, k, dt, beta, v_old)
% the pieces (owner K, source E, sign), their quadrature points x, the feet
% y = phi(x) and grad phi (see the header)
old = geo_old.mesh;
new = geo_new.mesh;
n_new = size(new.elements, 2);
[B, x3] = get_affine_map(new);
area_K = abs(B(1, 1, :).*B(2, 2, :) - B(1, 2, :).*B(2, 1, :))/2;
area_K = area_K(:).';
% feet of the vertices (reference vertices (1, 0), (0, 1), (0, 0))
rv = [1 0 0; 0 1 0];
Pv = dubiner_basis(k, rv);                                        % 3 x d2
d2 = size(Pv, 2);
X = reshape(new.coordinates(1:2, new.elements), 2, 3, n_new);
bx = Pv*beta(1:d2, :);
by = Pv*beta(d2 + 1:end, :);
b1 = reshape(permute(cat(3, bx, by), [3 1 2]), 2, []);           % 2 x 3 n_new
xv = reshape(X, 2, []);
y = xv - dt*b1;
for i = 1:60
    E = locate(geo_old, y);
    b0 = eval_vector(geo_old, k, v_old, E, y);
    y_new = xv - dt/2*(b1 + b0);
    change = max(abs(y_new - y), [], "all");
    y = y_new;
    if change < 1e-15, break, end
end
Yt = reshape(y, 2, 3, n_new);
% affine phi_K(x) = y_3 + A_K (x - x_3), A_K = [y1 - y3, y2 - y3] B^-1
Dy = [Yt(:, 1, :) - Yt(:, 3, :), Yt(:, 2, :) - Yt(:, 3, :)];
A = pagemrdivide(Dy, B);
det_A = reshape(A(1, 1, :).*A(2, 2, :) - A(1, 2, :).*A(2, 1, :), 1, []);
if any(det_A <= 0)
    error("sl_hyperelastic_solid_march:fold", "The affine feet fold %d element(s).", nnz(det_A <= 0));
end
Ai = pageinv(A);
T_old = reshape(old.coordinates(1:2, old.elements), 2, 3, []);
[K, Kp] = overlapping_boxes(Yt, T_old);
[P, n] = clip_polygons(counterclockwise(Yt(:, :, K)), counterclockwise(T_old(:, :, Kp)));
pair = zeros(1, 0);  a = zeros(2, 0);  b = a;  c = a;
for j = 2:size(P, 2) - 1
    has = find(n > j);
    pair = [pair, has]; %#ok<AGROW>
    a = [a, reshape(P(:, 1, has), 2, [])]; %#ok<AGROW>
    b = [b, reshape(P(:, j, has), 2, [])]; %#ok<AGROW>
    c = [c, reshape(P(:, j + 1, has), 2, [])]; %#ok<AGROW>
end
area_y = ((b(1, :) - a(1, :)).*(c(2, :) - a(2, :)) - (b(2, :) - a(2, :)).*(c(1, :) - a(1, :)))/2;
owner = K(pair);
keep = area_y > 1e-14*area_K(owner).*det_A(owner);
[pair, a, b, c, area_y, owner] = deal(pair(keep), a(:, keep), b(:, keep), c(:, keep), area_y(keep), owner(keep));
src = Kp(pair);
% back to K: x = x_3 + A^-1 (y - y_3)
back = @(z) reshape(pagemtimes(Ai(:, :, owner), reshape(z - reshape(Yt(:, 3, owner), 2, []), 2, 1, [])), ...
    2, []) + x3(:, owner);
[a, b, c] = deal(back(a), back(b), back(c));
area_x = area_y./det_A(owner);
covered = accumarray(owner(:), area_x(:), [n_new, 1]).';
uncovered = max(area_K - covered, 0)./area_K;
% extrapolation source: largest overlap, else the nearest old centroid
E_rem = zeros(1, n_new);
ov = accumarray([owner(:), src(:)], area_x(:), [n_new, size(T_old, 3)], [], [], true);
[m, i_max] = max(ov, [], 2);
E_rem(m > 0) = i_max(m > 0);
alone = find(E_rem == 0);
if ~isempty(alone)
    cK = reshape(mean(Yt(:, :, alone), 2), 2, []);
    E_rem(alone) = locate_nearest(geo_old, cK);
end
R = find(uncovered > 1e-12);
inR = ismember(owner, R);
[rs, w] = rule(k + 1);
nq = numel(w);
tri_points = @(a, b, c) reshape(a, 2, 1, []) + reshape(b - a, 2, 1, []).*rs(1, :) ...
    + reshape(c - a, 2, 1, []).*rs(2, :);
xK = pagemtimes(B(:, :, R), rs) + reshape(x3(:, R), 2, 1, []);
tr.x = cat(3, tri_points(a, b, c), xK, tri_points(a(:, inR), b(:, inR), c(:, inR)));
tr.owner = [owner, R, owner(inR)];
tr.src = [src, E_rem(R), E_rem(owner(inR))];
tr.W = w.*[area_x, area_K(R), -area_x(inR)];                       % nq x n_pieces
tr.uncovered = uncovered;
tr.geo_old = geo_old;
tr.geo_new = geo_new;
% true feet at the points, the source polynomial for v^n(y)
xs = reshape(tr.x, 2, []);
own = repelem(tr.owner, nq);
srs = repelem(tr.src, nq);
[b1, Gb1] = eval_vector(geo_new, k, beta, own, xs);
y = xs - dt*b1;
for i = 1:60
    b0 = eval_vector(geo_old, k, v_old, srs, y);
    y_new = xs - dt/2*(b1 + b0);
    change = max(abs(y_new - y), [], "all");
    y = y_new;
    if change < 1e-15, break, end
end
[~, Gb0] = eval_vector(geo_old, k, v_old, srs, y);
% grad phi = (I + dt/2 Gb0)^-1 (I - dt/2 Gb1), 2 x 2 x n pages
P1 = reshape(eye(2), 4, 1) + dt/2*Gb0;
P2 = reshape(eye(2), 4, 1) - dt/2*Gb1;
tr.J = reshape(pagemldivide(reshape(P1, 2, 2, []), reshape(P2, 2, 2, [])), 4, []);
tr.y = y;
tr.nq = nq;
end

function out = transport_state(tr, st, k)
% the transported fields of st (v, u, g in P_k^2 when present; L in S_h,
% weighted by grad phi)
go = tr.geo_old;
gn = tr.geo_new;
names = intersect(["v", "u", "g"], string(fieldnames(st)).', "stable");
Ud = cell2mat(cellfun(@(f) st.(f), cellstr(names), UniformOutput = false).');
d2 = nchoosek(k + 2, 2);
d3 = nchoosek(k + 3, 2);
n_old = size(go.mesh.elements, 2);
UL = reshape(pagemtimes(go.PtoD, reshape(st.L, [], 1, n_old)), [], n_old);   % 4*d3 x n_old
md = size(Ud, 1)/d2;
n_new = size(gn.mesh.elements, 2);
mom_d = zeros(d2*md, n_new);
mom_L = zeros(d3*4, n_new);
nq = tr.nq;
n_pieces = numel(tr.owner);
chunk = max(1, floor(4e5/nq));
for i0 = 1:chunk:n_pieces
    i = i0:min(i0 + chunk - 1, n_pieces);
    own = repelem(tr.owner(i), nq);
    srs = repelem(tr.src(i), nq);
    x = reshape(tr.x(:, :, i), 2, []);
    y = tr.y(:, (i0 - 1)*nq + 1:(i(end))*nq);
    J = tr.J(:, (i0 - 1)*nq + 1:(i(end))*nq);
    W = reshape(tr.W(:, i), 1, []);
    n = numel(own);
    inc = sparse(1:n, own, 1, n, n_new);
    rn = ref_coords(gn, own, x);
    ro = ref_coords(go, srs, y);
    Tk = dubiner_basis(k, rn).';                                     % d2 x n
    Sk = dubiner_basis(k, ro).';
    for c = 1:md
        f = sum(Sk.*Ud((c - 1)*d2 + (1:d2), srs), 1);
        mom_d((c - 1)*d2 + (1:d2), :) = mom_d((c - 1)*d2 + (1:d2), :) + (Tk.*(W.*f))*inc;
    end
    Tk1 = dubiner_basis(k + 1, rn).';
    Sk1 = dubiner_basis(k + 1, ro).';
    Lv = zeros(4, n);
    for c = 1:4
        Lv(c, :) = sum(Sk1.*UL((c - 1)*d3 + (1:d3), srs), 1);
    end
    F = mm2(Lv, J);
    for c = 1:4
        mom_L((c - 1)*d3 + (1:d3), :) = mom_L((c - 1)*d3 + (1:d3), :) + (Tk1.*(W.*F(c, :)))*inc;
    end
end
vol = gn.mesh.volume(:).';
Vd = kron(eye(md), gn.Mref)\mom_d./vol;
VL = kron(eye(4), gn.Mref1)\mom_L./vol;
out = struct();
for c = 1:numel(names)
    out.(names(c)) = Vd((c - 1)*2*d2 + (1:2*d2), :);
end
out.L = reshape(pagemtimes(gn.DtoP, reshape(VL, [], 1, n_new)), [], n_new);
end

function C = mm2(A, B)
% 2 x 2 products of column-major 4 x n arrays
C = [A(1, :).*B(1, :) + A(3, :).*B(2, :); A(2, :).*B(1, :) + A(4, :).*B(2, :); ...
     A(1, :).*B(3, :) + A(3, :).*B(4, :); A(2, :).*B(3, :) + A(4, :).*B(4, :)];
end

function r = ref_coords(geo, E, x)
r = reshape(pagemtimes(geo.Ginv(:, :, E), reshape(x - geo.x3(:, E), 2, 1, [])), 2, []);
end

function [V, GV] = eval_vector(geo, k, U, E, x)
% values (2 x n) and gradients (4 x n, column-major) of the P_k^2 field U
% at the points x, each with the polynomial of its element E
d2 = nchoosek(k + 2, 2);
r = ref_coords(geo, E, x);
if nargout > 1
    [P, Pr, Ps] = dubiner_basis(k, r);
else
    P = dubiner_basis(k, r);
end
V = [sum(P.'.*U(1:d2, E), 1); sum(P.'.*U(d2 + 1:end, E), 1)];
if nargout > 1
    G = geo.Ginv(:, :, E);
    GV = zeros(4, numel(E));
    for c = 1:2
        Uc = U((c - 1)*d2 + (1:d2), E);
        dr = sum(Pr.'.*Uc, 1);
        ds = sum(Ps.'.*Uc, 1);
        GV(c, :) = dr.*reshape(G(1, 1, :), 1, []) + ds.*reshape(G(2, 1, :), 1, []);       % d/dx
        GV(c + 2, :) = dr.*reshape(G(1, 2, :), 1, []) + ds.*reshape(G(2, 2, :), 1, []);   % d/dy
    end
end
end

function E = locate(geo, y)
% the element of geo's mesh containing each point (the nearest centroid
% outside the mesh)
E = pointLocation(geo.tri, y.');
out = isnan(E);
E(out) = locate_nearest(geo, y(:, out));
E = E(:).';
end

function E = locate_nearest(geo, y)
E = zeros(1, size(y, 2));
for i0 = 1:256:size(y, 2)
    i = i0:min(i0 + 255, size(y, 2));
    [~, E(i)] = min((geo.centroid(1, :).' - y(1, i)).^2 + (geo.centroid(2, :).' - y(2, i)).^2, [], 1);
end
end

function out = transport_check(geo0, geo1, k, dt, ex, options)
% one transport of the t = 0 patch level with the affine velocity
% v = V0 + omega (x - c0)^perp, against the exactly composed fields
om = options.transport_check;
c0 = options.center;
V0 = options.patch_velocity;
Wm = om*[0 -1; 1 0];
vel = @(x, y) [V0(1) - om*(y - c0(2)); V0(2) + om*(x - c0(1))];
st = initial_level(geo0, k, ex, 0);
st.v =project_vector(geo0.mesh, k, @(x, y) reshape(vel(x(:).', y(:).').', [], 1));
beta = project_vector(geo1.mesh, k, @(x, y) reshape(vel(x(:).', y(:).').', [], 1));
tr = transport_setup(geo0, geo1, k, dt, beta, st.v);
Xt = transport_state(tr, st, k);
% exact feet: (I + dt/2 W) y = (I - dt/2 W) x - dt (V0 - W c0)
Mi = inv(eye(2) + dt/2*Wm);
Jp = Mi*(eye(2) - dt/2*Wm);
b0 = -dt*Mi*(V0 - Wm*c0);
at = @(h, x, y) h(Jp(1, 1)*x + Jp(1, 2)*y + b0(1), Jp(2, 1)*x + Jp(2, 2)*y + b0(2), 0);
vx = project_vector(geo1.mesh, k, @(x, y) at(@(a, b, t) reshape(vel(a.', b.').', [], 1), x, y));
ux = project_vector(geo1.mesh, k, @(x, y) at(ex.u, x, y));
Lx = project_peers(geo1, k, @(x, y) reshape(mm2(reshape(at(ex.grad_u, x, y), [], 4).', ...
    repmat(Jp(:), 1, numel(x))).', [], 1));
rel = [l2_dubiner(geo1, k, Xt.v - vx)/l2_dubiner(geo1, k, vx), ...
    l2_peers(geo1, Xt.L - Lx)/l2_peers(geo1, Lx), l2_dubiner(geo1, k, Xt.u - ux)/l2_dubiner(geo1, k, ux)];
out = struct("transport_error", max(rel), "transport_errors", rel, "uncovered", max(tr.uncovered), ...
    "n_pieces", numel(tr.owner));
end

function e = l2_dubiner(geo, k, U)
d2 = nchoosek(k + 2, 2);
e = 0;
for c = 0:size(U, 1)/d2 - 1
    X = U(c*d2 + (1:d2), :);
    e = e + sum(sum(X.*(geo.Mref*X), 1).*geo.mesh.volume(:).');
end
e = sqrt(e);
end

function e = l2_peers(geo, S)
n_el = size(S, 2);
e = sqrt(sum(pagemtimes(reshape(S, [], 1, n_el), "transpose", ...
    pagemtimes(geo.MP, reshape(S, [], 1, n_el)), "none"), "all"));
end

%% ------------------------------------------------------- geometry, data
function geo = geometry(mesh, k, options)
% the operators of one carve
geo.mesh = mesh;
n_el = size(mesh.elements, 2);
d1 = k + 1;
d2 = nchoosek(k + 2, 2);
d3 = nchoosek(k + 3, 2);
lin = get_matrices(LinearElasticityProblem(k, young_modulus = options.E, ...
    poisson_ratio = options.nu), mesh, k);
p = size(lin{1}, 1) - 3*d2;
geo.p = p;
geo.tS = lin{1}(1:p, 1:p, :);                       % (A sigma, tau)
geo.tU = lin{1}(1:p, p + (1:2*d2), :);              % (u, div tau)
geo.tG = lin{1}(1:p, p + 2*d2 + (1:d2), :);         % (gamma, tau)
geo.tH = lin{2}(1:p, :, :);                         % <u_hat, tau n>
geo.MP = PEERS_matrix(mesh, k);
geo.M = mass_matrix(mesh, k);
[geo.tPP, nxLP, nyLP, geo.tLP, geo.tLL] = boundary_matrix(mesh, k);
O = zeros(3*d1, d2, n_el);
o = zeros(3*d1, p - 4*d2, n_el);
geo.Sn = [nxLP, O, nyLP, O, o; O, nxLP, O, nyLP, o];
n_dof = d1*size(mesh.faces, 2);
face_dof = reshape(((mesh.boundaries(:) - 1)*d1 + (1:d1)).', 3*d1, n_el);
geo.idx = [face_dof; face_dof + n_dof];
geo.n_dof = n_dof;
% volume rule of the nonlinear terms, PEERS values there
[rs, w] = rule(ceil(3*(k + 1)/2) + 1);
geo.Pp = peers_points(mesh, k, rs);
geo.wv = w.*mesh.volume(:).';
geo.rs = rs;
% P_k+1 Dubiner <-> S_h on every element (S_h is inside P_k+1^(2x2))
[rq, wq] = rule(k + 2);
Pq = peers_points(mesh, k, rq);
D1 = dubiner_basis(k + 1, rq);
nq = numel(wq);
geo.Mref = dubiner_basis(k, rq).'*(wq.*dubiner_basis(k, rq));
geo.Mref1 = D1.'*(wq.*D1);
C = zeros(p, 4*d3, n_el);
for c = 1:4
    C(:, (c - 1)*d3 + (1:d3), :) = pagemtimes(Pq((c - 1)*nq + (1:nq), :, :), "transpose", ...
        (wq.*D1).*reshape(mesh.volume, 1, 1, []), "none");
end
geo.DtoP = pagemldivide(geo.MP, C);
geo.PtoD = pagemtimes(kron(eye(4), inv(geo.Mref1)), pagetranspose(C))./reshape(mesh.volume, 1, 1, []);
[geo.Ginv, geo.x3] = get_affine_map(mesh, inverse = true);
geo.tri = triangulation(double(mesh.elements(1:3, :).'), mesh.coordinates(1:2, :).');
geo.centroid = reshape(mean(reshape(mesh.coordinates(1:2, mesh.elements), 2, 3, []), 2), 2, []);
end

function st = initial_level(geo, k, ex, t)
st.v = project_vector(geo.mesh, k, @(x, y) ex.v(x, y, t));
st.L = project_peers(geo, k, @(x, y) ex.grad_u(x, y, t));
st.u = project_vector(geo.mesh, k, @(x, y) ex.u(x, y, t));
st.g = project_vector(geo.mesh, k, @(x, y) ex.g(x, y, t));
end

function data = step_data(geo, k, ex, t)
% sources and boundary data at t
mesh = geo.mesh;
data.f = moments(mesh, k, @(x, y) ex.f(x, y, t));
data.g = project_vector(mesh, k, @(x, y) ex.g(x, y, t));
data.gradg = peers_moments(geo, @(x, y) ex.grad_g(x, y, t));
sp = LinearElasticityProblem(k);
sp.Data.boundary_conditions = boundary_condition("Dirichlet", ...
    x_data = @(x, y) row(ex.v(x(:), y(:), t), 1, size(x)), ...
    y_data = @(x, y) row(ex.v(x(:), y(:), t), 2, size(x)));
[data.lam, data.dir] = get_boundary_conditions(sp, mesh);
data.dir = data.dir(:);
end

function r = row(s, i, sz)
n = numel(s)/2;
r = reshape(s((i - 1)*n + (1:n)), sz);
end

function F = peers_moments(geo, h)
% (h, tau)_K for the S_h basis, h 2 x 2 valued (column-major blocks)
[xs, ys] = points(geo.mesh, geo.rs);
s = unstack(h(xs(:), ys(:)), xs);                                  % nq*4 x n_el
F = pagemtimes(geo.Pp, "transpose", reshape(repmat(geo.wv, 4, 1).*s, [], 1, size(xs, 2)), "none");
end

function L = project_peers(geo, k, h) %#ok<INUSD>
n_el = size(geo.mesh.elements, 2);
L = reshape(pagemldivide(geo.MP, peers_moments(geo, h)), [], n_el);
end

function B = blk2(A)
Z = zeros(size(A, 1), size(A, 2), size(A, 3));
B = [A, Z; Z, A];
end

%% ---------------------------------------------- carve, fields, measures
function grid = local_background(N, center, radius)
% the cells of StructuredSimplicialGrid(N, N, "crossed") on [0, 1]^2 that
% lie within two cells of the disc's bounding box
lo = max(floor((center - radius)*N) - 2, 0);
hi = min(ceil((center + radius)*N) + 2, N);
grid = StructuredSimplicialGrid(hi(1) - lo(1), hi(2) - lo(2), lo(1)/N, hi(1)/N, lo(2)/N, hi(2)/N, ...
    style = "crossed");
end

function mesh = carve_solid(background, curve, snap)
% the inside of the circle, cut elements discarded, interface vertices
% snapped onto it; every boundary face Dirichlet; tau = 1 (the face
% matrices are scaled by tau_v in the rows)
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

function ex = exact_fields(options, amplitude, lambda, mu)
% handles of (x, y, t): v, u, f, g, dtv = D_t v (2 x n), sigma, grad u,
% grad g, grad_vg = grad v + grad g (column-major 2 x 2)
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
    case "quadratic"
        u = 0.03*[sin(2*x + y); cos(x - 2*y)] + t*0.02*[cos(3*y); sin(2*x)] ...
            + t^2*0.02*[x*y; x - y^2];
        v = 0.2*[cos(x + 2*y); sin(2*x - y)] + t*0.1*[sin(y); cos(x)] + t^2*0.1*[x*y; x^2];
    case "patch"
        V0 = sym(options.patch_velocity);
        Xp = xx - V0*t;
        U = sym([0; 0]);
        i = 0;
        for a = 0:options.patch_degree
            for b = 0:options.patch_degree - a
                i = i + 1;
                U = U + Xp(1)^a*Xp(2)^b*[cos(3 + i); sin(5 + 2*i)];
            end
        end
        u = 0.05*U*(1 + t^2);
        v = V0 + 0*x;
end
G = jacobian(u, xx);
g = simplify(diff(u, t) + G*v - v);
if options.fields == "standard", g = sym([0; 0]); end
if options.law == "linear"
    sigma = mu*(G + G.') + lambda*trace(G)*eye(2);
else
    a = eye(2) - G;
    j = det(a);
    sigma = j*(mu*(inv(a.'*a) - eye(2)) - lambda*log(j)*eye(2));
end
dtv = diff(v, t) + jacobian(v, xx)*v;
f = options.rho*dtv - [diff(sigma(1, 1), x) + diff(sigma(1, 2), y); diff(sigma(2, 1), x) + diff(sigma(2, 2), y)];
Gg = jacobian(g, xx);
Gv = jacobian(v, xx);
ex.v = stacked(v, x, y, t);
ex.u = stacked(u, x, y, t);
ex.f = stacked(f, x, y, t);
ex.g = stacked(g, x, y, t);
ex.dtv = stacked(dtv, x, y, t);
ex.sigma = stacked(sigma(:), x, y, t);
ex.grad_u = stacked(G(:), x, y, t);
ex.grad_g = stacked(Gg(:), x, y, t);
ex.grad_vg = stacked(Gv(:) + Gg(:), x, y, t);
end

function h = stacked(e, x, y, t)
hs = arrayfun(@(ei) matlabFunction(ei, Vars = [x, y, t], Optimize = false), e, UniformOutput = false);
h = @(xx, yy, tt) cell2mat(cellfun(@(f) f(xx, yy, tt) + 0*xx, hs, UniformOutput = false));
end

function [rs, w] = rule(q)
[rs, w] = get_2d_reference_quadrature(q);
w = w(:)/sum(w);
end

function [xs, ys] = points(mesh, rs)
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
n_c = numel(s)/numel(xs);
V = zeros(n_c*size(D, 2), size(xs, 2));
for i = 1:n_c
    V((i - 1)*size(D, 2) + (1:size(D, 2)), :) = fit*reshape(s((i - 1)*numel(xs) + (1:numel(xs))), nq, []);
end
end

function F = moments(mesh, k, h)
[rs, w] = rule(ceil(3*k/2) + 1);
[xs, ys] = points(mesh, rs);
s = h(xs(:), ys(:));
D = dubiner_basis(k, rs);
n = numel(xs);
F = [(w.*D).'*reshape(s(1:n), numel(w), []); (w.*D).'*reshape(s(n + 1:2*n), numel(w), [])];
F = reshape(F.*mesh.volume(:).', [], 1, size(xs, 2));
end

function Pp = peers_points(mesh, k, rs)
% PEERS basis at the reference points rs, (nq*4) x nP x n_el, rows q
% fastest then the components 11, 21, 12, 22
G = get_affine_map(mesh, inverse = true);
P = PEERS_basis(k, rs(1, :), rs(2, :), G);
[nq, nP] = size(P, 1:2);
Pp = reshape(permute(P, [1 3 4 2 5]), nq*4, nP, []);
end

function V = unstack(s, xs)
[nq, n] = size(xs);
V = reshape(permute(reshape(s, nq, n, []), [1 3 2]), [], n);
end

function E = energy(geo, k, prm, V, L)
% 1/2 rho |v|^2 + 1/2 (C eps(L), eps(L))
d2 = nchoosek(k + 2, 2);
nq = size(geo.Pp, 1)/4;
Lq = reshape(pagemtimes(geo.Pp, reshape(L, [], 1, size(L, 2))), nq*4, []);
h = @(i) Lq((i - 1)*nq + (1:nq), :);
s12 = (h(2) + h(3))/2;
tr = h(1) + h(4);
es = sum(geo.wv.*(prm.mu*(h(1).^2 + h(4).^2 + 2*s12.^2) + prm.lambda/2*tr.^2), "all");
ev = sum(sum(V(1:d2, :).*(geo.Mref*V(1:d2, :)) + V(d2 + 1:end, :).*(geo.Mref*V(d2 + 1:end, :)), 1) ...
    .*geo.mesh.volume(:).');
E = prm.rho*ev/2 + es;
end

function [e, nrm] = errors(geo, k, st, S, ex, t)
% L2 errors and exact norms [v, L, sigma, u_P] at time t
mesh = geo.mesh;
d2 = nchoosek(k + 2, 2);
[rs, w] = rule(k + 3);
[xs, ys] = points(mesh, rs);
D = dubiner_basis(k, rs);
W = w.*mesh.volume(:).';
nq = numel(w);
Pp = peers_points(mesh, k, rs);
pv = @(A) reshape(pagemtimes(Pp, reshape(A, [], 1, size(A, 2))), nq*4, []);
sq = @(a, nc) reshape(sum(reshape(a, nq, nc, []).^2, 2), nq, []);
num = {kron(eye(2), D)*st.v, pv(st.L), pv(S), kron(eye(2), D)*st.u};
exact = {unstack(ex.v(xs(:), ys(:), t), xs), unstack(ex.grad_u(xs(:), ys(:), t), xs), ...
    unstack(ex.sigma(xs(:), ys(:), t), xs), unstack(ex.u(xs(:), ys(:), t), xs)};
nc = [2, 4, 4, 2];
[e, nrm] = deal(zeros(1, 4));
for i = 1:4
    e(i) = sqrt(sum(W.*sq(num{i} - exact{i}, nc(i)), "all"));
    nrm(i) = sqrt(sum(W.*sq(exact{i}, nc(i)), "all"));
end
end

function b = field_bounds(mesh, k, ex, t, c_s)
[rs, ~] = rule(k + 3);
[xs, ys] = points(mesh, rs);
n = numel(xs);
G = reshape(ex.grad_u(xs(:), ys(:), t), n, 4);
v = reshape(ex.v(xs(:), ys(:), t), n, 2);
b = [max(sqrt(sum(G.^2, 2))), max(sqrt(sum(v.^2, 2)))/c_s, ...
    min((1 - G(:, 1)).*(1 - G(:, 4)) - G(:, 2).*G(:, 3))];
end

%% ----------------------------------------------------- polygon clipping
function T = counterclockwise(T)
signed = (T(1, 2, :) - T(1, 1, :)).*(T(2, 3, :) - T(2, 1, :)) ...
       - (T(2, 2, :) - T(2, 1, :)).*(T(1, 3, :) - T(1, 1, :));
flip = signed < 0;
T(:, [2 3], flip) = T(:, [3 2], flip);
end

function [P, n] = clip_polygons(subject, clip)
% Sutherland-Hodgman, vectorized over pages: subject ∩ clip, both
% counterclockwise triangles (2 x 3 x n); n the vertex counts
n_pages = size(subject, 3);
MAX = 6;
P = zeros(2, MAX, n_pages);
P(:, 1:3, :) = subject;
n = 3*ones(1, n_pages);
pages = 1:n_pages;
for e = 1:3
    MAX_in = size(P, 2);
    MAX = max(MAX_in, 2*max([n, 0]));
    a = reshape(clip(:, e, :), 2, []);
    b = reshape(clip(:, mod(e, 3) + 1, :), 2, []);
    side = @(p) (b(1, :) - a(1, :)).*(p(2, :) - a(2, :)) - (b(2, :) - a(2, :)).*(p(1, :) - a(1, :));
    Q = zeros(2, MAX, n_pages);
    n_new = zeros(1, n_pages);
    for i = 1:max([n, 0])
        active = i <= n;
        j_prev = (i - 1)*ones(1, n_pages);
        if i == 1, j_prev = max(n, 1); end
        current = reshape(P(:, i, :), 2, []);
        prior = P(:, (pages - 1)*MAX_in + j_prev);
        d_cur = side(current);
        d_pri = side(prior);
        in_cur = d_cur >= 0;
        in_pri = d_pri >= 0;
        crossing = active & (in_cur ~= in_pri);
        lam = d_pri./(d_pri - d_cur + ~crossing);
        X = prior + (current - prior).*lam;
        n_new = n_new + crossing;
        Q = put(Q, X, n_new, crossing);
        keep = active & in_cur;
        n_new = n_new + keep;
        Q = put(Q, current, n_new, keep);
    end
    P = Q;
    n = n_new;
end
end

function Q = put(Q, X, slot, mask)
[~, MAX, n_pages] = size(Q);
p = find(mask);
Q = reshape(Q, 2, MAX*n_pages);
Q(:, (p - 1)*MAX + slot(p)) = X(:, p);
Q = reshape(Q, 2, MAX, n_pages);
end

function [K, Kp] = overlapping_boxes(T, S)
% the pairs (K, K') of triangles T(:, :, K), S(:, :, K') whose bounding
% boxes meet (uniform bins of the largest box side)
lo = reshape(min(T, [], 2), 2, []);   hi = reshape(max(T, [], 2), 2, []);
lo_s = reshape(min(S, [], 2), 2, []); hi_s = reshape(max(S, [], 2), 2, []);
width = max([hi - lo, hi_s - lo_s], [], "all");
slack = 1e-12*width;
origin = min(lo, [], 2) - slack;
top = max(hi, [], 2) + slack;
n_bins = max(ceil((top - origin)/width), 1);
near = find(all(hi_s + slack >= origin & lo_s - slack <= top, 1));
A = bin_incidence(lo - slack, hi + slack, origin, width, n_bins);
A_s = bin_incidence(lo_s(:, near) - slack, hi_s(:, near) + slack, origin, width, n_bins);
[K, j] = find(A*A_s.');
Kp = near(j);
K = reshape(K, 1, []);
Kp = reshape(Kp, 1, []);
meets = all(lo(:, K) <= hi_s(:, Kp) + slack & lo_s(:, Kp) <= hi(:, K) + slack, 1);
K = K(meets);
Kp = Kp(meets);
end

function A = bin_incidence(lo, hi, origin, width, n_bins)
bin = @(x) min(max(floor((x - origin)/width) + 1, 1), n_bins);
i0 = bin(lo);
i1 = bin(hi);
rows = [];  cols = [];
for dx = 0:max(i1(1, :) - i0(1, :))
    for dy = 0:max(i1(2, :) - i0(2, :))
        has = find(i0(1, :) + dx <= i1(1, :) & i0(2, :) + dy <= i1(2, :));
        rows = [rows, has]; %#ok<AGROW>
        cols = [cols, i0(1, has) + dx + (i0(2, has) + dy - 1)*n_bins(1)]; %#ok<AGROW>
    end
end
A = sparse(rows, cols, 1, size(lo, 2), prod(n_bins));
end
