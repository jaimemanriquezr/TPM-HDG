function out = velocity_stress_solid_march(k, N, amplitude, dt, n_steps, options)
%VELOCITY_STRESS_SOLID_MARCH  Velocity-stress elastic solid on a moving snapped carve.
%
%   out = velocity_stress_solid_march(k, N, A, dt, n_steps, rotation = 0, ...
%       scheme = "trapezoid", fields = "standard", homogeneous = false, ...)
%
%   The solid alone, Eulerian: the disc of the given radius moves rigidly,
%   phi_t(X) = c(t) + Q(theta(t)) (X - c0), c(t) = c0 + A sin(2 pi t) e_x,
%   theta(t) = rotation sin(2 pi t), Q the rotation matrix. At every t_n
%   the solid mesh is the inside of the circle Sigma(t_n) carved from the
%   background StructuredSimplicialGrid(N, N, "crossed") on [0, 1]^2 (cut
%   elements discarded; only the cells near the disc are built) with its interface vertices snapped onto Sigma(t_n)
%   (snap_to_curve, normal_sign = -1); every boundary face is Dirichlet.
%   Unknowns: the velocity v (relative to the rigid motion) and the
%   stress sigma, small strain, with the co-rotational material derivative
%   D v = d_t v + beta.grad v - W v, D sigma = d_t sigma + beta.grad sigma
%   - W sigma + sigma W, beta the rigid velocity field, W = theta' J:
%     rho D v - div sigma = f,   A D sigma - eps(v) = 0,
%   A = C^-1 the compliance. The state of a step is (sigma, v) only.
%   One step t_n -> t_n+1, R = Q(theta_n+1 - theta_n), foot(x) = c_n +
%   R' (x - c_n+1):
%   1. transport to the new carve: v~ = Pi_V[R v_n(foot)], sigma~ the
%      A-orthogonal projection onto the weakly symmetric PEERS stresses of
%      R sigma_n(foot) R' (find sigma~ in S_h, xi in P_k skew:
%      (A sigma~, tau) + (xi, tau) = (A R sigma_n(foot) R', tau),
%      (sigma~, eta) = 0). The old fields are read through
%      LinearElasticityProblem.translated_field_between on the old carve
%      moved rigidly onto the new frame: exact overlaps, the polynomial of
%      one old element extrapolated where foot(x) leaves the old carve.
%      The stress goes through its exact P_k+1 Dubiner representation.
%   2. scheme = "trapezoid": the PEERS HDG solve for the midpoint values
%      (sigma', v', gamma', v^') on the new carve,
%        (2/dt)(A sigma', tau) + (v', div tau) + (gamma', tau) - <v^', tau n> = (2/dt)(A sigma~, tau)
%        (2 rho/dt)(v', w) - (div sigma', w) + <S(v' - v^'), w> = (2 rho/dt)(v~, w) + (f', w)
%        (sigma', eta) = 0,  v^' = g' on the boundary,
%      f' = (R f(foot, t_n) + f(., t_n+1))/2 and g' likewise from the
%      exact velocity; then sigma_n+1 = 2 sigma' - sigma~, v_n+1 = 2 v' - v~
%      (Crank-Nicolson). gamma' is the midpoint skew multiplier and is
%      not carried. scheme = "bdf2": the same solve with 3/(2 dt) in place
%      of 2/dt, history (4 X~_n - X~_n-1)/3 (X~_n-1 the previous step's
%      transported level, transported once more) and data at t_n+1; the
%      first step is a trapezoid step.
%   With c = 2/dt (or 3/(2 dt)) the solve is a LinearElasticityProblem
%   with Young modulus E/c (compliance c A), density c rho, Newmark with
%   beta = 1 and time_step = 1 and the history v~ as displacement (the
%   reaction c rho (v, w) and the right-hand side c rho (v~, w)), and
%   Data.constitutive_source = c (A sigma~, tau).
%   fields = "standard": body-frame displacement U(X, t) = u_S(phi_t(X), t),
%   u_S = options.u_solid; v = Q d_t U, sigma = Q C eps_X(U) Q' at X =
%   phi_t^-1(x) (with rotation = 0, v = D_t u_S and sigma = C eps(u_S)).
%   fields = "patch": body-frame V(X) of degree k and Sigma(X, t) = S0(X)
%   + t C eps_X(V), S0 symmetric of degree k: reproduced to round-off.
%   homogeneous = true: f = 0, g = 0 and a random initial (sigma, v) of
%   unit energy (seed).
%   initial_perturbation = [seed, amplitude]: adds a random (sigma, v) of
%   energy amplitude^2 times the initial energy.
%
%   out: e_sigma, e_v (L2 errors at the final time, NaN when homogeneous),
%   n_sigma, n_v (L2 norms of the exact fields there); energy (n_steps+1
%   x 1, 1/2 (rho |v|^2 + (A sigma, sigma)) after every solve, the first
%   the initial one), energy_transported (n_steps x 1, after each
%   transport); uncovered (n_steps x 2: the largest fraction of a new
%   element and the total area not covered by the moved old carve);
%   step_errors (n_steps x 2, record_errors only); the final sigma, v,
%   mesh; seconds.
%
%See also: LinearElasticityProblem.translated_field_between, snap_to_curve
arguments
    k (1,1) double
    N (1,1) double
    amplitude (1,1) double
    dt (1,1) double
    n_steps (1,1) double
    options.rotation (1,1) double = 0
    options.scheme (1,1) string {mustBeMember(options.scheme, ["trapezoid", "bdf2"])} = "trapezoid"
    options.fields (1,1) string {mustBeMember(options.fields, ["standard", "patch"])} = "standard"
    options.u_solid (2,1) string = ["0.1*sin(2*x + y)*cos(t)"; "0.1*cos(x - y)*(1 + t^2)"]
    options.homogeneous (1,1) logical = false
    options.seed (1,1) double = 1
    options.initial_perturbation double = []
    options.record_errors (1,1) logical = false
    options.center (2,1) double = [0.5; 0.5]
    options.radius (1,1) double = 0.2
    options.snap (1,1) logical = true
    % continue a trapezoid march from step start_step with the state
    % start_state = struct(sigma, v) (out.sigma, out.v of the run to step
    % start_step); the records then cover the steps after it
    options.start_step (1,1) double {mustBeInteger, mustBeNonnegative} = 0
    options.start_state struct = struct.empty
    options.rho (1,1) double = 1
    options.E (1,1) double = 10
    options.nu (1,1) double = 0.3
end
wall = tic;
rho = options.rho;
[lambda, mu] = LinearElasticityProblem.get_lame_parameters(options.E, options.nu);
c1 = 1/(2*mu);
c2 = lambda/(4*mu*(lambda + mu));
c0 = options.center;
cen = @(t) c0 + [amplitude*sin(2*pi*t); 0];
ang = @(t) options.rotation*sin(2*pi*t);
rot = @(a) [cos(a), -sin(a); sin(a), cos(a)];
moving = amplitude ~= 0 || options.rotation ~= 0;
carve = @(t) carve_solid(local_background(N, cen(t), options.radius), ...
    CircularCurve(cen(t), options.radius), options.snap);
ex = exact_fields(k, options, amplitude, lambda, mu);
d2 = nchoosek(k + 2, 2);
sp0 = LinearElasticityProblem(k);   % only the object translated_field_between needs

n0 = options.start_step;
assert(n0 == 0 || (options.scheme == "trapezoid" && ~isempty(options.start_state)), ...
    "velocity_stress_solid_march:start", "A restart needs the trapezoid scheme and a start_state.");
mesh = carve(n0*dt);
ops = stress_operators(mesh, k, c1, c2);
if n0 > 0
    [S, V] = deal(options.start_state.sigma, options.start_state.v);
elseif options.homogeneous
    stream = RandStream("mt19937ar", Seed = options.seed);
    [S, V] = random_state(stream, mesh, k, ops, rho);
else
    V = project_vector(mesh, k, @(x, y) ex.v(x, y, 0));
    S = project_stress(mesh, k, ops, c1, c2, @(x, y) ex.sigma(x, y, 0));
end
if ~isempty(options.initial_perturbation)
    stream = RandStream("mt19937ar", Seed = options.initial_perturbation(1));
    [dS, dV] = random_state(stream, mesh, k, ops, rho);
    a = options.initial_perturbation(2)*sqrt(energy(mesh, k, ops, rho, S, V));
    S = S + a*dS;
    V = V + a*dV;
end
E = nan(n_steps - n0 + 1, 1);
E_t = nan(n_steps - n0, 1);
E(1) = energy(mesh, k, ops, rho, S, V);
uncovered = zeros(n_steps - n0, 2);
step_errors = nan(n_steps - n0, 2);
previous = [];   % bdf2: the last transported level on the current carve
for n = n0 + 1:n_steps
    t0 = (n - 1)*dt;
    t1 = n*dt;
    R = rot(ang(t1) - ang(t0));
    [c_old, c_new] = deal(cen(t0), cen(t1));
    if moving
        mesh_new = carve(t1);
        ops = stress_operators(mesh_new, k, c1, c2);
        [Vt, St, uncovered(n - n0, :)] = transport(sp0, mesh, mesh_new, k, V, S, c_old, c_new, R, ops, c1, c2);
        if ~isempty(previous)
            [previous.V, previous.S] = transport(sp0, mesh, mesh_new, k, previous.V, previous.S, ...
                c_old, c_new, R, ops, c1, c2);
        end
        mesh = mesh_new;
    else
        [Vt, St] = deal(V, S);
    end
    E_t(n - n0) = energy(mesh, k, ops, rho, St, Vt);
    bdf2 = options.scheme == "bdf2" && ~isempty(previous);
    if bdf2
        c = 3/(2*dt);
        [Zv, Zs] = deal((4*Vt - previous.V)/3, (4*St - previous.S)/3);
        f = at_time(ex.f, t1);
        g = at_time(ex.v, t1);
    else
        c = 2/dt;
        [Zv, Zs] = deal(Vt, St);
        f = averaged(ex.f, t0, t1, R, c_old, c_new);
        g = averaged(ex.v, t0, t1, R, c_old, c_new);
    end
    if options.homogeneous
        zero = @(x, y) 0*x;
        [f, g] = deal({zero; zero});
    end
    sp = LinearElasticityProblem(k, density = c*rho, time_step = 1, time_scheme = "newmark", ...
        newmark_beta = 1, newmark_gamma = 0.5, young_modulus = options.E/c, ...
        poisson_ratio = options.nu, source_term = f);
    sp.Data.boundary_conditions = boundary_condition("Dirichlet", x_data = g{1}, y_data = g{2});
    sp.Data.state = struct("solid", struct("displacement", Zv, "velocity", 0, "acceleration", 0));
    sp.Data.constitutive_source = c*reshape(pagemtimes(ops.T, reshape(Zs, [], 1, size(Zs, 2))), ...
        size(Zs));
    sol = solve(sp, mesh, k);
    Ys = dofs_of(sol.Unknowns.field(1).dof);
    Yv = dofs_of(sol.Unknowns.field(2).dof);
    if options.scheme == "bdf2"
        previous = struct("V", Vt, "S", St);
    end
    if bdf2
        [S, V] = deal(Ys, Yv);
    else
        [S, V] = deal(2*Ys - St, 2*Yv - Vt);
    end
    E(n - n0 + 1) = energy(mesh, k, ops, rho, S, V);
    if options.record_errors && ~options.homogeneous
        [step_errors(n - n0, 1), step_errors(n - n0, 2)] = errors(mesh, k, S, V, ex, t1);
    end
end
T = n_steps*dt;
[e_sigma, e_v, n_sigma, n_v] = deal(NaN);
if ~options.homogeneous
    [e_sigma, e_v, n_sigma, n_v] = errors(mesh, k, S, V, ex, T);
end
out = struct("e_sigma", e_sigma, "e_v", e_v, "n_sigma", n_sigma, "n_v", n_v, "energy", E, ...
    "energy_transported", E_t, "uncovered", uncovered, "step_errors", step_errors, ...
    "sigma", S, "v", V, "mesh", mesh, "n_elements", size(mesh.elements, 2), ...
    "seconds", toc(wall));
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
% snapped onto it; every boundary face Dirichlet
if snap
    [~, mesh] = background.split(curve, max_normal_deviation = 90);
    mesh = snap_to_curve(mesh, curve, normal_sign = -1);
else
    [~, mesh] = background.split(curve);
end
is_boundary = accumarray(mesh.boundaries(:), 1, [size(mesh.faces, 2), 1]) == 1;
mesh.faces(3, is_boundary) = mesh.boundary_tags("Dirichlet");
end

function ex = exact_fields(k, options, amplitude, lambda, mu)
% handles v(x, y, t) {2}, sigma(x, y, t) {2 x 2}, f(x, y, t) {2}
syms x y t X Y real
c0 = options.center;
c = c0 + [amplitude*sin(2*pi*t); 0];
th = options.rotation*sin(2*pi*t);
Q = [cos(th), -sin(th); sin(th), cos(th)];
C_eps = @(U) 2*mu*sym_eps(U, X, Y) + lambda*trace(sym_eps(U, X, Y))*eye(2);
switch options.fields
    case "standard"
        u = str2sym(options.u_solid(:));
        phi = c + Q*([X; Y] - c0);
        U = subs(u, [x; y], phi);
        Vb = diff(U, t);
        Sb = C_eps(U);
    case "patch"
        Vb = sym(zeros(2, 1));
        S0 = sym(zeros(2, 2));
        i = 0;
        for a = 0:k
            for b = 0:k - a
                i = i + 1;
                m = X^a*Y^b;
                Vb = Vb + m*[cos(1 + i); sin(2 + 3*i)];
                s = [cos(3 + i), sin(5 + 2*i); 0, cos(7 + 5*i)];
                S0 = S0 + m*(s + triu(s, 1).');
            end
        end
        Sb = S0 + t*C_eps(Vb);
end
Fb = options.rho*diff(Vb, t) - [diff(Sb(1, 1), X) + diff(Sb(1, 2), Y); ...
    diff(Sb(2, 1), X) + diff(Sb(2, 2), Y)];
Xx = c0 + Q.'*([x; y] - c);
to_x = @(e) subs(e, [X; Y], Xx);
v = to_x(Q*Vb);
sigma = to_x(Q*Sb*Q.');
f = to_x(Q*Fb);
handle = @(e) vectorized(matlabFunction(e, Vars = [x, y, t]));
hv = arrayfun(handle, v, UniformOutput = false);
hs = arrayfun(handle, sigma, UniformOutput = false);
hf = arrayfun(handle, f, UniformOutput = false);
ex.v = @(xx, yy, tt) cellfun(@(h) h(xx, yy, tt), hv, UniformOutput = false);
ex.sigma = @(xx, yy, tt) cellfun(@(h) h(xx, yy, tt), hs, UniformOutput = false);
ex.f = @(xx, yy, tt) cellfun(@(h) h(xx, yy, tt), hf, UniformOutput = false);
end

function E = sym_eps(U, X, Y)
G = [diff(U(1), X), diff(U(1), Y); diff(U(2), X), diff(U(2), Y)];
E = (G + G.')/2;
end

function h = vectorized(f)
h = @(x, y, t) f(x, y, t) + 0*x;
end

function h = at_time(F, t)
% {2} handles (x, y) of a vector field at time t
h = {@(x, y) pick(F(x, y, t), 1); @(x, y) pick(F(x, y, t), 2)};
end

function h = averaged(F, t0, t1, R, c_old, c_new)
% {2} handles (x, y) of (R F(foot(x), t0) + F(x, t1))/2
h = cell(2, 1);
for i = 1:2
    h{i} = @(x, y) (rotated_at_foot(F, x, y, t0, R, c_old, c_new, i) + pick(F(x, y, t1), i))/2;
end
end

function w = rotated_at_foot(F, x, y, t, R, c_old, c_new, i)
fx = c_old(1) + R(1, 1)*(x - c_new(1)) + R(2, 1)*(y - c_new(2));
fy = c_old(2) + R(1, 2)*(x - c_new(1)) + R(2, 2)*(y - c_new(2));
Ff = F(fx, fy, t);
w = R(i, 1)*Ff{1} + R(i, 2)*Ff{2};
end

function v = pick(c, i)
v = c{i};
end

function D = dofs_of(D)
if iscell(D), D = cat(1, D{:}); end
end

function [rs, w] = rule(q)
[rs, w] = get_2d_reference_quadrature(q);
w = w(:)/sum(w);
end

function ops = stress_operators(mesh, k, c1, c2)
% per element: T = (A sigma, tau) matrix, the weak symmetry block MA and
% the saddle matrix of the A-orthogonal weakly symmetric projection
[M_P, M_tr, M_A] = PEERS_matrix(mesh, k);
ops.T = c1*M_P - c2*M_tr;
nP = size(M_P, 1);
d2 = size(M_A, 2);
ops.K = [ops.T, M_A; permute(M_A, [2 1 3]), zeros(d2, d2, size(M_A, 3))];
ops.nP = nP;
end

function S = sym_projection(ops, rhs)
% sigma in S_h with (A sigma, tau) + (xi, tau) = rhs(tau), (sigma, eta) = 0
n_el = size(rhs, 2);
d2 = size(ops.K, 1) - ops.nP;
X = pagemldivide(ops.K, cat(1, reshape(rhs, ops.nP, 1, n_el), zeros(d2, 1, n_el)));
S = reshape(X(1:ops.nP, 1, :), ops.nP, n_el);
end

function Pp = peers_points(mesh, k, rs)
% PEERS basis at the reference points rs of every element, as
% (nq*4) x nP x n_el, rows q fastest then the components (1,1), (2,1),
% (1,2), (2,2)
G = get_affine_map(mesh, inverse = true);
P = PEERS_basis(k, rs(1, :), rs(2, :), G);
[nq, nP] = size(P, 1:2);
Pp = reshape(permute(P, [1 3 4 2 5]), nq*4, nP, []);
end

function vals = stress_values(mesh, k, S, rs)
% (nq*4) x n_el values of the PEERS fields S at rs
Pp = peers_points(mesh, k, rs);
vals = reshape(pagemtimes(Pp, reshape(S, size(S, 1), 1, [])), size(Pp, 1), []);
end

function rhs = stress_moments(mesh, k, Avals, rs, w)
% (A sigma, tau_j)_K from the values of A sigma ((nq*4) x n_el) at rs
Pp = peers_points(mesh, k, rs);
W = repmat(w, 4, 1);
rhs = reshape(pagemtimes(Pp, "transpose", reshape(W.*Avals, [], 1, size(Avals, 2)), "none"), ...
    size(Pp, 2), []).*mesh.volume;
end

function Avals = compliance(vals, c1, c2)
% A sigma at the points, values (nq*4) x n_el in the component order
nq = size(vals, 1)/4;
tr = vals(1:nq, :) + vals(3*nq + (1:nq), :);
Avals = c1*vals;
Avals(1:nq, :) = Avals(1:nq, :) - c2*tr;
Avals(3*nq + (1:nq), :) = Avals(3*nq + (1:nq), :) - c2*tr;
end

function vals = rotate_tensor(vals, R)
% R sigma R' at every point
nq = size(vals, 1)/4;
s = @(i) vals((i - 1)*nq + (1:nq), :);
S = {s(1), s(3); s(2), s(4)};   % S{a, b} = sigma_ab
out = cell(2, 2);
for a = 1:2
    for b = 1:2
        out{a, b} = 0;
        for p = 1:2
            for q = 1:2
                out{a, b} = out{a, b} + R(a, p)*S{p, q}*R(b, q);
            end
        end
    end
end
vals = [out{1, 1}; out{2, 1}; out{1, 2}; out{2, 2}];
end

function S = project_stress(mesh, k, ops, c1, c2, sigma)
% A-orthogonal weakly symmetric projection of the handle sigma(x, y) {2 x 2}
[rs, w] = rule(k + 3);
xy = mesh.get_2d_physical(repmat(rs, 1, 1, size(mesh.elements, 2)), 1:size(mesh.elements, 2));
s = sigma(reshape(xy(1, :, :), numel(w), []), reshape(xy(2, :, :), numel(w), []));
vals = [s{1, 1}; s{2, 1}; s{1, 2}; s{2, 2}];
S = sym_projection(ops, stress_moments(mesh, k, compliance(vals, c1, c2), rs, w));
end

function V = project_vector(mesh, k, v)
% L2 projection onto P_k of the handle v(x, y) {2}
[rs, w] = rule(k + 3);
xy = mesh.get_2d_physical(repmat(rs, 1, 1, size(mesh.elements, 2)), 1:size(mesh.elements, 2));
s = v(reshape(xy(1, :, :), numel(w), []), reshape(xy(2, :, :), numel(w), []));
D = dubiner_basis(k, rs);
fit = (D.'*(w.*D)) \ (D.'.*w.');
V = [fit*s{1}; fit*s{2}];
end

function [S, V] = random_state(stream, mesh, k, ops, rho)
% random (sigma, v) of unit energy, sigma weakly symmetric
n_el = size(mesh.elements, 2);
d2 = nchoosek(k + 2, 2);
V = randn(stream, 2*d2, n_el);
S0 = randn(stream, ops.nP, n_el);
S = sym_projection(ops, reshape(pagemtimes(ops.T, reshape(S0, ops.nP, 1, n_el)), ops.nP, n_el));
a = 1/sqrt(energy(mesh, k, ops, rho, S, V));
S = a*S;
V = a*V;
end

function E = energy(mesh, k, ops, rho, S, V)
% 1/2 (rho |v|^2 + (A sigma, sigma))
[rs, w] = rule(k + 1);
D = dubiner_basis(k, rs);
M = D.'*(w.*D);
d2 = size(M, 1);
ev = 0;
for c = 0:1
    Vc = V(c*d2 + (1:d2), :);
    ev = ev + sum(sum(Vc.*(M*Vc), 1).*mesh.volume);
end
n_el = size(S, 2);
es = sum(pagemtimes(reshape(S, [], 1, n_el), "transpose", ...
    pagemtimes(ops.T, reshape(S, [], 1, n_el)), "none"), "all");
E = (rho*ev + es)/2;
end

function [Vt, St, uncovered] = transport(sp0, old, new, k, V, S, c_old, c_new, R, ops, c1, c2)
% (sigma, v) of old carried by the rigid step onto new (see the header)
moved = old;
X = old.coordinates;
moved.coordinates = c_new + R*(X - c_old);
nn = reshape(old.normals, 2, []);
moved.normals = reshape(R*nn, size(old.normals));
% velocity: P_k Dubiner, then rotated components
[Vt, info] = sp0.translated_field_between(moved, new, k, k, V, [0; 0]);
d2 = nchoosek(k + 2, 2);
v1 = Vt(1:d2, :);
v2 = Vt(d2 + (1:d2), :);
Vt = [R(1, 1)*v1 + R(1, 2)*v2; R(2, 1)*v1 + R(2, 2)*v2];
area = new.volume(:);
uncovered = [max(info.uncovered), sum(info.uncovered(:).*area)];
% stress: exact P_k+1 Dubiner representation on old, carried, rotated,
% A-orthogonal weakly symmetric projection on new
[rs, w] = rule(k + 2);
Dk1 = dubiner_basis(k + 1, rs);
fit = (Dk1.'*(w.*Dk1)) \ (Dk1.'.*w.');
vals = stress_values(old, k, S, rs);
nq = numel(w);
F = zeros(4*size(fit, 1), size(S, 2));
d3 = size(fit, 1);
for i = 1:4
    F((i - 1)*d3 + (1:d3), :) = fit*vals((i - 1)*nq + (1:nq), :);
end
assert(norm(kron(eye(4), Dk1)*F - vals, "fro") <= 1e-9*max(1, norm(vals, "fro")), ...
    "velocity_stress_solid_march:peers", "PEERS stresses not in P_k+1.");
Ft = sp0.translated_field_between(moved, new, k + 1, k + 1, F, [0; 0]);
vals_t = kron(eye(4), Dk1)*Ft;
vals_t = rotate_tensor(vals_t, R);
St = sym_projection(ops, stress_moments(new, k, compliance(vals_t, c1, c2), rs, w));
end

function [e_s, e_v, n_s, n_v] = errors(mesh, k, S, V, ex, t)
% L2 errors and exact norms of sigma and v at time t
[rs, w] = rule(k + 3);
n_el = size(mesh.elements, 2);
xy = mesh.get_2d_physical(repmat(rs, 1, 1, n_el), 1:n_el);
xs = reshape(xy(1, :, :), numel(w), []);
ys = reshape(xy(2, :, :), numel(w), []);
s = ex.sigma(xs, ys, t);
s = [s{1, 1}; s{2, 1}; s{1, 2}; s{2, 2}];
sh = stress_values(mesh, k, S, rs);
W4 = repmat(w, 4, 1);
e_s = sqrt(sum(sum(W4.*(sh - s).^2, 1).*mesh.volume));
n_s = sqrt(sum(sum(W4.*s.^2, 1).*mesh.volume));
v = ex.v(xs, ys, t);
D = dubiner_basis(k, rs);
d2 = size(D, 2);
vh = {D*V(1:d2, :), D*V(d2 + (1:d2), :)};
e_v = sqrt(sum(sum(w.*((vh{1} - v{1}).^2 + (vh{2} - v{2}).^2), 1).*mesh.volume));
n_v = sqrt(sum(sum(w.*(v{1}.^2 + v{2}.^2), 1).*mesh.volume));
end
