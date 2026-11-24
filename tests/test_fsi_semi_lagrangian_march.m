% test_fsi_semi_lagrangian_march.m
% eulerian_solid_march with transport = "foot": the solid Newmark history
% carried along the feet of each step (transport_solid_state,
% foot_projection) for a reference moving rigidly, a disc translating and
% rotating about its centre.
%   1. transport_solid_state with the exact rigid foot (rotation by 0.3 and
%      translation by 0.05 e_x, so that elements switch on and off) on the
%      snapped circle of radius 0.2, crossed grid N = 16, k = 1, 2: a history
%      whose every level is a global polynomial p of degree k (volume,
%      traces, the lifted stress and rotation) goes to the history of p o
%      foot on the new carve (1e-10, every field). The same with grad foot
%      replaced by I is reported: the stress and rotation levels are then
%      not those of p o foot. The transported momentum flux P~ (the fourth
%      output) in place of the stress gives the consistent acceleration
%      (LinearElasticityProblem.consistent_acceleration, zero source) of the
%      momentum identity at the feet, the projection of
%      (div C sym grad p)(foot(x))/rho (1e-10 relative to the largest
%      expected dof, or absolute below 1); the same with the transported
%      stress level is reported (it is not, for k = 2).
%   2. uniform field: (a) the foot transport with foot(x) = x - shift
%      equals translate_state + remap_solid_state on global polynomials
%      (1e-10): flat y = 0.45 -> 0.55 on the "left" grid N = 8 and the
%      circle shifted by 0.05 e_x on the crossed grid N = 16, snapped,
%      k = 1, 2; (b) the march with the transport omitted is the march with
%      transport = "translate" (identical dofs, 3 steps); (c) on a discrete
%      state (three steps of the smooth march on the moving circle, k = 1,
%      N = 16) the largest dof difference between the two transports per
%      field class, over the largest dof of the class, separately on the
%      kept elements and on the moved ones (where the rules differ: E2
%      extrapolates the source's translate, the foot evaluates at the
%      feet). With foot_integration = "quadrature" at quadrature_degree
%      2k + 2 and 2k + 6 (a quadrature error of a discontinuous integrand
%      on the kept elements), bounded loosely (below 1). With "clipped"
%      (exact, as E2's clipped translate): the volume levels on the kept
%      elements agree to 1e-10; the stress and rotation differ there by
%      construction (the foot transports C sym(G) and skew(G) of the
%      gradient G = c1 sigma - c2 tr(sigma) I + gamma, E2 sigma and gamma:
%      they differ by the skew part of the weakly symmetric PEERS stress),
%      and agree to 1e-10 on the same state with each stress replaced by
%      its symmetric Dubiner part (bubbles dropped).
%   3. the march is exact on a solution compatible with the rigid motion:
%      circle, A = 0.1, Theta = 0.5, snapped, foot_rule = "exact", crossed
%      N = 8, dt = 0.05, 6 steps, u_S = U0(X) + t U1(X) + t^2 U2(X), X =
%      phi_t^-1(x), U_j of degree k (P_k in x at every t, quadratic in t
%      along the characteristics), stream function x y (1 + t), pressure
%      (x + y)(1 + t): worst L2 error per step below 1e-9, k = 1, 2, and the
%      interface rebuilt from u* and from u_h on the rotated reference
%      (measure_interface) at the exact material points (1e-9). The same
%      bounds for the march with consistent_flux = "transported" (the
%      consistent acceleration from the transported momentum flux, kept
%      under rotation) and skeleton_history = "ustar_displacement" (the
%      displacement trace from u* at the feet), with history_velocity =
%      "ustar_twostep" and history_acceleration = "consistent_kept", its
%      interface rebuilt on the material pre-image of Sigma(t_n)
%      (interface_reference = "preimage"), where the distance to the exact
%      set is bounded as well (1e-9). The
%      default march with foot_rule = "midpoint", foot_velocity =
%      "extrapolated" is reported.
%   4. smooth data over a period: circle, k = 1, N = 8, A = 0.1, Theta = 0.5,
%      dt = 0.05, T = 1, default foot options: bounded (max|u_h| < 10), and
%      the worst error over the steps of (sigma_F, u_F, sigma_S, u_S) below
%      [1.3 1.25 5.5 1.3] times that of the march at Theta = 0 with
%      transport = "translate" (each with its default history, see
%      eulerian_solid_march's transport): 1.25 to 1.45 times the measured
%      ratios [1 1 4.3 0.906]. The foot march at Theta = 0 and the march with
%      foot_rule = "exact" are reported.
clc; clear;
warning("off", "MATLAB:nearlySingularMatrix");   % trace extension on sliver patches

n_fail = 0;
classes = struct("volume", ["displacement", "velocity", "acceleration"], ...
    "skeleton", ["displacement_skeleton", "velocity_skeleton", "acceleration_skeleton"], ...
    "stress", ["stress", "stress_velocity", "stress_acceleration"], ...
    "rotation", ["rotation", "rotation_velocity", "rotation_acceleration"]);
identity = @(x) repmat(eye(2), 1, 1, size(x, 2));

%% 1. foot transport exact on P_k, rigid rotation and translation
alpha = 0.3;  c_old = [0.5; 0.5];  c_new = [0.55; 0.5];
R = [cos(alpha), sin(alpha); -sin(alpha), cos(alpha)];         % Q(-alpha)
foot = @(x) c_old + R*(x - c_new);
grad_foot = @(x) repmat(R, 1, 1, size(x, 2));
background = StructuredSimplicialGrid(16, 16, "crossed");
for k = 1:2
    [sp, p, q] = polynomial_data(k);
    [old, old_keys, old_faces] = carve(background, CircularCurve(c_old, 0.2));
    [new, new_keys, new_faces] = carve(background, CircularCurve(c_new, 0.2));
    s = levels(sp, old, k, p, q);
    expected = levels(sp, new, k, composed(p, foot), composed(q, foot));
    [t, switched, ~, flux] = transport_solid_state(sp, s, old, new, k, old_keys, new_keys, old_faces, ...
        new_faces, foot, grad_foot);
    d = class_differences(t, expected, classes);
    t_I = transport_solid_state(sp, s, old, new, k, old_keys, new_keys, old_faces, new_faces, ...
        foot, identity);
    d_I = class_differences(t_I, expected, classes);
    fprintf("1. rigid foot k = %d: max difference volume %.2e skeleton %.2e stress %.2e rotation %.2e;" + ...
        " switched on/off %d/%d\n", k, d, switched);
    fprintf("   grad foot = I: volume %.2e skeleton %.2e stress %.2e rotation %.2e\n", d_I);
    n_fail = n_fail + check(max(d) <= 1e-10 && all(switched > 0), sprintf("rigid foot transport k = %d", k));
    % the consistent acceleration from the transported momentum flux
    f = stress_divergence(sp, k);
    expected_a = L2proj2D(new, {@(x, y) f(1) + 0*x; @(x, y) f(2) + 0*x}, k)/sp.Data.density;
    flux_state = struct("solid", t);
    a_level = consistent_acceleration(sp, new, k, flux_state);
    flux_state.solid.stress = flux;
    a_flux = consistent_acceleration(sp, new, k, flux_state);
    largest = max(abs(expected_a), [], "all");
    e_flux = max(abs(a_flux - expected_a), [], "all");
    fprintf("   consistent acceleration at the feet: from the momentum flux %.2e, from the stress level %.2e" + ...
        " (largest expected dof %.2e)\n", e_flux, max(abs(a_level - expected_a), [], "all"), largest);
    n_fail = n_fail + check(e_flux <= 1e-10*max(1, largest), sprintf("transported momentum flux k = %d", k));
end

%% 2a. uniform foot against translate_state + remap_solid_state
cases = {StructuredSimplicialGrid(8, 8, "left"), {FlatCurveHorizontal(0.45), FlatCurveHorizontal(0.55)}, ...
    [0; 0.1], "flat"; background, {CircularCurve([0.5; 0.5], 0.2), CircularCurve([0.55; 0.5], 0.2)}, ...
    [0.05; 0], "circle"};
for c = 1:size(cases, 1)
    [grid, curves, shift, label] = cases{c, :};
    for k = 1:2
        [sp, p, q] = polynomial_data(k);
        [old, old_keys, old_faces] = carve(grid, curves{1});
        [new, new_keys, new_faces] = carve(grid, curves{2});
        s = levels(sp, old, k, p, q);
        e2 = e2_transport(sp, s, old, new, k, old_keys, new_keys, old_faces, new_faces, shift);
        t = transport_solid_state(sp, s, old, new, k, old_keys, new_keys, old_faces, new_faces, ...
            @(x) x - shift, identity);
        d = class_differences(t, e2, classes);
        fprintf("2a. %s k = %d, uniform foot vs translate + remap: volume %.2e skeleton %.2e stress %.2e" + ...
            " rotation %.2e\n", label, k, d);
        n_fail = n_fail + check(max(d) <= 1e-10, sprintf("uniform foot = E2 on P_k, %s k = %d", label, k));
    end
end

%% 2b. the default transport is translate
uu = "sin(pi*x)*sin(pi*y)*cos(t)";  pp = "cos(pi*x)*(y - 0.5)*(1 + t)";
us = ["0.1*sin(2*x + y)*cos(t)"; "0.1*cos(x - y)*(1 + t^2)"];
flat_grid = StructuredSimplicialGrid(8, 8, "left");
out_default = eulerian_solid_march(1, flat_grid, 0.1, 0.05, 3, uu, pp, us);
out_translate = eulerian_solid_march(1, flat_grid, 0.1, 0.05, 3, uu, pp, us, transport = "translate");
d = max(abs([dofs(out_default.fs) - dofs(out_translate.fs); dofs(out_default.ss) - dofs(out_translate.ss)]));
fprintf("2b. default vs transport = ""translate"": max dof difference %.2e\n", d);
n_fail = n_fail + check(d == 0, "the default transport is not the translate");

%% 2c. a discrete state: foot quadrature against the exact clipped translate
k = 1;  A = 0.1;  dt = 0.05;  n_steps = 3;
[~, out] = evalc("eulerian_solid_march(k, background, A, dt, n_steps, uu, pp, us, geometry = ""circle"")");
sp = LinearElasticityProblem(k, young_modulus = 10, poisson_ratio = 0.3, time_scheme = "newmark", ...
    density = 1);
shift = A*(sin(2*pi*(n_steps + 1)*dt) - sin(2*pi*n_steps*dt))*[1; 0];
[new, new_keys, new_faces] = carve(background, CircularCurve([0.5 + A*sin(2*pi*(n_steps + 1)*dt); 0.5], 0.2));
old = out.solid_mesh;
[e2, ~, moved] = e2_transport(sp, out.state.solid, old, new, k, out.keys.solid, new_keys, ...
    out.keys.solid_faces, new_faces, shift);
kept = setdiff(1:size(new.elements, 2), moved);
zero = structfun(@(f) 0*f, e2, UniformOutput = false);
scale = class_differences(e2, zero, classes);
fprintf("2c. discrete state, k = 1, shift %.4f, %d kept and %d moved elements; max |dof| volume %.2e" + ...
    " skeleton %.2e stress %.2e rotation %.2e\n", shift(1), numel(kept), numel(moved), scale);
for qd = [2*k + 2, 2*k + 6]
    t = transport_solid_state(sp, out.state.solid, old, new, k, out.keys.solid, new_keys, ...
        out.keys.solid_faces, new_faces, @(x) x - shift, identity, quadrature_degree = qd, ...
        foot_integration = "quadrature");
    d_kept = class_differences(t, e2, classes, kept);
    d_moved = class_differences(t, e2, classes, moved);
    fprintf("    quadrature_degree %d, max dof difference / max |dof|: kept elements volume %.2e stress %.2e" + ...
        " rotation %.2e; moved elements volume %.2e stress %.2e rotation %.2e; faces skeleton %.2e\n", qd, ...
        d_kept([1 3 4])./scale([1 3 4]), d_moved([1 3 4])./scale([1 3 4]), d_kept(2)/scale(2));
    n_fail = n_fail + check(max([d_kept, d_moved]./[scale, scale]) < 1, sprintf("discrete state, quadrature_degree %d", qd));
end
t = transport_solid_state(sp, out.state.solid, old, new, k, out.keys.solid, new_keys, ...
    out.keys.solid_faces, new_faces, @(x) x - shift, identity);
d_kept = class_differences(t, e2, classes, kept);
d_moved = class_differences(t, e2, classes, moved);
fprintf("    clipped, max dof difference / max |dof|: kept elements volume %.2e stress %.2e" + ...
    " rotation %.2e; moved elements volume %.2e stress %.2e rotation %.2e; faces skeleton %.2e\n", ...
    d_kept([1 3 4])./scale([1 3 4]), d_moved([1 3 4])./scale([1 3 4]), d_kept(2)/scale(2));
n_fail = n_fail + check(d_kept(1)/scale(1) <= 1e-10 && max([d_kept, d_moved]./[scale, scale]) < 1, ...
    "discrete state, clipped: volume on the kept elements");
s_sym = symmetric_stress(out.state.solid, k, classes.stress);
e2 = e2_transport(sp, s_sym, old, new, k, out.keys.solid, new_keys, out.keys.solid_faces, new_faces, shift);
t = transport_solid_state(sp, s_sym, old, new, k, out.keys.solid, new_keys, ...
    out.keys.solid_faces, new_faces, @(x) x - shift, identity);
d_kept = class_differences(t, e2, classes, kept);
d_moved = class_differences(t, e2, classes, moved);
fprintf("    clipped, symmetric Dubiner stresses: kept elements volume %.2e stress %.2e rotation %.2e;" + ...
    " moved elements volume %.2e stress %.2e rotation %.2e\n", ...
    d_kept([1 3 4])./scale([1 3 4]), d_moved([1 3 4])./scale([1 3 4]));
n_fail = n_fail + check(max(d_kept([1 3 4])./scale([1 3 4])) <= 1e-10, ...
    "discrete state, clipped, symmetric stresses: kept elements");

%% 3. exact solution compatible with the rigid motion
A = 0.1;  Theta = 0.5;  dt = 0.05;  n_steps = 6;
grid8 = StructuredSimplicialGrid(8, 8, "crossed");
for k = 1:2
    u_poly = rigid_polynomial(k, A, Theta, [0.5; 0.5]);
    march = @(varargin) eulerian_solid_march(k, grid8, A, dt, n_steps, "x*y*(1 + t)", "(x + y)*(1 + t)", ...
        u_poly, "geometry", "circle", "transport", "foot", "rotation", Theta, "verbose", true, varargin{:});
    timer = tic;
    [~, out] = evalc("march(""foot_rule"", ""exact"", ""measure_interface"", true)");
    wall = toc(timer);
    fprintf("3. polynomial, rigid motion, k = %d, exact foot: worst L2 error per step %s; switched" + ...
        " (F on/off, S on/off) %s; %.1f s\n", k, mat2str(out.errors.', 3), ...
        mat2str(sum(out.history(:, 4:7), 1)), wall);
    point_error = max(out.interface(:, [3 9]), [], 1);
    fprintf("   interface rebuilt on the rotated reference: worst point error from u* %.2e, from u_h %.2e\n", ...
        point_error);
    n_fail = n_fail + check(max(out.errors) <= 1e-9, sprintf("polynomial march exact, k = %d", k));
    n_fail = n_fail + check(max(point_error) <= 1e-9, sprintf("polynomial march k = %d: interface points", k));
    n_fail = n_fail + check(all(sum(out.history(:, [6 7]), 1) > 0), ...
        sprintf("polynomial march k = %d: solid elements do not switch on and off", k));
    flux_options = {"foot_rule", "exact", "measure_interface", true, "history_velocity", "ustar_twostep", ...
        "history_acceleration", "consistent_kept", "consistent_flux", "transported", ...
        "skeleton_history", "ustar_displacement", "interface_reference", "preimage"};
    [~, out] = evalc("march(flux_options{:})");
    point_error = max(out.interface(:, [3 9]), [], 1);
    distance = max(out.interface(:, [2 8]), [], 1);
    fprintf("   transported flux, u* displacement trace: worst L2 error per step %s; interface on the" + ...
        " pre-image: points from u* %.2e, from u_h %.2e, distance from u* %.2e, from u_h %.2e\n", ...
        mat2str(out.errors.', 3), point_error, distance);
    n_fail = n_fail + check(max(out.errors) <= 1e-9, sprintf("polynomial march, transported flux, k = %d", k));
    n_fail = n_fail + check(max([point_error, distance]) <= 1e-9, ...
        sprintf("polynomial march, transported flux, k = %d: interface on the pre-image", k));
    [~, out] = evalc("march(""foot_rule"", ""midpoint"", ""foot_velocity"", ""extrapolated"")");
    fprintf("   midpoint, extrapolated: worst L2 error per step %s\n", mat2str(out.errors.', 3));
end

%% 4. smooth data over a period
k = 1;  dt = 0.05;  n_steps = 20;
grid8 = StructuredSimplicialGrid(8, 8, "crossed");
timer = tic;
[~, out] = evalc("eulerian_solid_march(k, grid8, 0.1, dt, n_steps, uu, pp, us, geometry = ""circle"", " + ...
    "transport = ""foot"", rotation = 0.5, verbose = true)");
wall = toc(timer);
[~, ref] = evalc("eulerian_solid_march(k, grid8, 0.1, dt, n_steps, uu, pp, us, geometry = ""circle"", " + ...
    "verbose = true)");
[~, foot0] = evalc("eulerian_solid_march(k, grid8, 0.1, dt, n_steps, uu, pp, us, geometry = ""circle"", " + ...
    "transport = ""foot"", verbose = true)");
[~, exact] = evalc("eulerian_solid_march(k, grid8, 0.1, dt, n_steps, uu, pp, us, geometry = ""circle"", " + ...
    "transport = ""foot"", rotation = 0.5, foot_rule = ""exact"", verbose = true)");
worst = max(out.step_errors, [], 1);
ratio = worst ./ max(ref.step_errors, [], 1);
bound = [1.3, 1.25, 5.5, 1.3];
fprintf("4. smooth, circle, A = 0.1, Theta = 0.5, T = 1, worst over steps: |sF-s| %.3e |uF-u| %.3e" + ...
    " |sS-s| %.3e |uS-u| %.3e; %.1f s\n", worst, wall);
fprintf("   Theta = 0, translate: |sF-s| %.3e |uF-u| %.3e |sS-s| %.3e |uS-u| %.3e\n", ...
    max(ref.step_errors, [], 1));
fprintf("   Theta = 0, foot: |sF-s| %.3e |uF-u| %.3e |sS-s| %.3e |uS-u| %.3e\n", ...
    max(foot0.step_errors, [], 1));
fprintf("   Theta = 0.5, exact foot: |sF-s| %.3e |uF-u| %.3e |sS-s| %.3e |uS-u| %.3e\n", ...
    max(exact.step_errors, [], 1));
fprintf("   ratio %s (bounds %s)\n", mat2str(ratio, 3), mat2str(bound));
fprintf("   max |u_h| (fluid, solid) %.3g %.3g; switched (F on/off, S on/off) %s\n", ...
    max(out.history(:, 2:3), [], 1), mat2str(sum(out.history(:, 4:7), 1)));
n_fail = n_fail + check(all(isfinite([out.eF, out.eS])) && max(out.history(:, 2:3), [], "all") < 10, ...
    "smooth: the march is unbounded");
n_fail = n_fail + check(all(ratio < bound), ...
    "smooth: worst error over the steps of (sF, uF, sS, uS) above its bound times the translate march's");

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_fsi_semi_lagrangian_march:failed", "%d check(s) failed", n_fail);
end

function [sp, p, q] = polynomial_data(k)
% the solid and two vector polynomials of degree k (the E2 test's)
sp = LinearElasticityProblem(k, young_modulus = 10, poisson_ratio = 0.3, time_scheme = "newmark", ...
    density = 1);
p = {@(x, y) x.^2 - 2*x.*y + y; @(x, y) y.^2 + x.*y - x};
if k == 1, p = {@(x, y) x - 2*y + 1; @(x, y) 3*y + x}; end
q = {@(x, y) 2*x - y; @(x, y) x + 0.5*y - 1};
end

function f = stress_divergence(sp, k)
% div C sym grad p of polynomial_data's p (constant): mu lap p + (lambda +
% mu) grad div p
[lambda, mu] = LinearElasticityProblem.get_lame_parameters(sp.Data.young_modulus, sp.Data.poisson_ratio);
f = [0; 0];
if k == 2, f = mu*[2; 2] + (lambda + mu)*[3; 0]; end
end

function s = symmetric_stress(s, k, names)
% each PEERS stress level replaced by the symmetric part of its Dubiner
% modes (blocks E11, E21, E12, E22 of d2 rows), the bubbles dropped
d2 = nchoosek(k + 2, 2);
for name = names
    S = s.(name);
    off = (S(d2 + (1:d2), :) + S(2*d2 + (1:d2), :))/2;
    S(d2 + (1:d2), :) = off;
    S(2*d2 + (1:d2), :) = off;
    S(4*d2 + 1:end, :) = 0;
    s.(name) = S;
end
end

function g = composed(f, foot)
% the handles x -> f(foot(x))
g = cell(size(f));
for c = 1:numel(f)
    h = f{c};
    g{c} = @(x, y) at_foot(h, foot, x, y);
end
end

function v = at_foot(h, foot, x, y)
z = foot([x(:).'; y(:).']);
v = reshape(h(z(1, :), z(2, :)), size(x));
end

function s = levels(solid, mesh, k, p, q)
% displacement and acceleration from p, velocity from q
[u, u_hat, sigma, gamma] = project(solid, mesh, k, p);
[v, v_hat, sigma_v, gamma_v] = project(solid, mesh, k, q);
s = struct("displacement", u, "velocity", v, "acceleration", u, ...
    "displacement_skeleton", u_hat, "velocity_skeleton", v_hat, "acceleration_skeleton", u_hat, ...
    "stress", sigma, "rotation", gamma, "stress_velocity", sigma_v, "rotation_velocity", gamma_v, ...
    "stress_acceleration", sigma, "rotation_acceleration", gamma);
end

function [u, u_hat, sigma, gamma] = project(solid, mesh, k, f)
u = L2proj2D(mesh, f, k);
u_hat = L2projSkeleton1d(mesh, f, k);
[sigma, gamma] = solid.lift_stress(mesh, k, u, u_hat);
end

function [s, switched, moved] = e2_transport(sp, s, old, new, k, old_keys, new_keys, old_faces, new_faces, shift)
% translate_state + remap_solid_state, as the translate march applies them
state.solid = s;
state = translate_state(sp, old, state, shift, k);
[s, switched, moved] = remap_solid_state(sp, state.solid, old, new, k, old_keys, new_keys, old_faces, ...
    new_faces, untranslated = s, shift = shift);
end

function d = class_differences(a, b, classes, elements)
% the largest dof difference per field class (volume, skeleton, stress,
% rotation); with elements, the element fields on those columns only (the
% skeleton on every face)
d = zeros(1, 4);
names = string(fieldnames(classes)).';
for i = 1:4
    for n = classes.(names(i))
        D = abs(as_matrix(a.(n)) - as_matrix(b.(n)));
        if nargin > 3 && names(i) ~= "skeleton", D = D(:, elements); end
        d(i) = max([d(i); D(:)]);
    end
end
end

function A = as_matrix(A)
if iscell(A), A = cat(1, A{:}); end
end

function u = rigid_polynomial(k, A, Theta, c0)
% u_S(x, t) = U0(X) + t U1(X) + t^2 U2(X), X = phi_t^-1(x) = c0 + Q(-theta)(x
% - c(t)), c(t) = c0 + A sin(2 pi t) e_x, theta = Theta sin(2 pi t), U_j of
% degree k: strings in x, y, t
th = sprintf("(%.17g*sin(2*pi*t))", Theta);
cx = sprintf("(%.17g + %.17g*sin(2*pi*t))", c0(1), A);
cy = sprintf("%.17g", c0(2));
X1 = sprintf("(%.17g + cos(%s)*(x - %s) + sin(%s)*(y - %s))", c0(1), th, cx, th, cy);
X2 = sprintf("(%.17g - sin(%s)*(x - %s) + cos(%s)*(y - %s))", c0(2), th, cx, th, cy);
if k == 1
    U0 = ["0.1*" + X1 + " - 0.05*" + X2; "0.02 + 0.1*" + X2];
    U1 = ["0.05*" + X2; "-0.04*" + X1 + " + 0.01"];
    U2 = ["0.03*" + X1 + " + 0.02*" + X2; "-0.02*" + X2];
else
    U0 = ["0.1*" + X1 + "*" + X2 + " - 0.05*" + X2; "0.1*(" + X1 + "^2 - " + X2 + ")"];
    U1 = ["0.05*" + X2 + "^2"; "0.05*" + X1 + "*" + X2 + " + 0.02"];
    U2 = ["0.03*" + X1 + "^2 - 0.01*" + X2; "-0.02*" + X2 + "^2 + 0.01*" + X1];
end
u = "(" + U0 + ") + t*(" + U1 + ") + t^2*(" + U2 + ")";
end

function [mesh, keys, face_keys] = carve(background, curve)
% the snapped solid side and its unsnapped element and face keys, as
% eulerian_solid_march carves it
[~, mesh] = background.split(curve, max_normal_deviation = 90);
keys = squeeze(mean(reshape(mesh.coordinates(:, mesh.elements), 2, 3, []), 2)).';
face_keys = ((mesh.coordinates(:, mesh.faces(1, :)) + mesh.coordinates(:, mesh.faces(2, :)))/2).';
mesh = snap_to_curve(mesh, curve, normal_sign = -1);
end

function x = dofs(solution)
% all field and hybrid dofs of a solution, stacked
c = [{solution.Unknowns.field.dof}, {solution.Unknowns.hybrid.dof}];
c = cellfun(@(a) reshape(double(cat(1, a{:})), [], 1), cellfun(@cellify, c, UniformOutput = false), ...
    UniformOutput = false);
x = cat(1, c{:});
end

function c = cellify(x)
if iscell(x), c = x(:); else, c = {x}; end
end
