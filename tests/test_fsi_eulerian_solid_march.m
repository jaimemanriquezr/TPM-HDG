% test_fsi_eulerian_solid_march.m
% eulerian_solid_march (tests/helpers): fluid and solid both carved from
% one background grid at every step, Sigma(t) prescribed, the solid
% Newmark history translated along the characteristics of the known
% velocity d'(t) (translate_state) and remapped across carves
% (remap_solid_state), its volume levels replaced by default by the exact
% projection of the translated field (direct_translate_state).
%   1. remap_solid_state after translate_state is exact on polynomials of
%      degree k: the volume, skeleton, stress and rotation fields equal the
%      projections of the translated polynomial on the new carve (1e-10).
%      Flat Sigma from y = 0.45 to 0.55 on the "left" grid N = 8 (one row
%      of 2N solid elements switches on), k = 1, 2, unsnapped and snapped
%      (snapped also with the untranslated faces carried along the shift);
%      and the circle of radius 0.2 on the crossed grid N = 16, snapped
%      (snap_to_curve, normal_sign = -1), shifted by 0.03 e_x (elements
%      switch off only) and by 0.05 e_x (elements also switch on, with the
%      nearest old element as source on a curved carve), both with the
%      untranslated state given, so that the boundary velocity and
%      acceleration traces are slid along the shift.
%   2. amplitude 0 reproduces a march on one fixed carve (no remap, no
%      translation) dof for dof (1e-12): k = 1, N = 8, snapped, BDF2, with
%      the toolbox history (history_translate = skeleton_history =
%      "library", the translated acceleration, trapezoidal Newmark); the
%      default skeleton_history = "ustar_displacement" replaces the
%      displacement trace history and differs from it by 1.97e-2.
%   3. the march is exact on a solution every ingredient reproduces: flat,
%      y0 = 0.45, A = 0.1, N = 8, k = 2, dt = 0.05, 6 steps (Sigma crosses
%      the grid line y = 0.5: a row switches on in the solid and off in the
%      fluid), u_S = [x (y - d - 0.45); x^2 - (y - d)] (1 + t + t^2/2),
%      stream function x y (1 + t), pressure (x + y)(1 + t). Worst L2 error
%      per step below 1e-9, snapped and unsnapped.
%   4. smooth data, k = 1, N = 8, A = 0.1, dt = 0.05, 9 steps (T = 0.45:
%      d and d' nonzero, Sigma has crossed y = 0.5 up and back), snapped:
%      the march is bounded (max|u_h| < 10), elements switch on and off on
%      both sides, and the worst L2 error over the steps of each of sigma_F,
%      u_F, sigma_S, u_S is below [4 1.65 2.5 1.15] times that of the same
%      march at A = 0 (one fixed carve): about 1.25 times the measured
%      ratios [3.17 1.3 2 0.913]; with history_velocity = "newmark" they
%      are [3.17 1.31 2 1.01].
%   5. the snapped circle of radius 0.2 at rest (A = 0) with time-independent
%      data, crossed grid N = 32, k = 1, dt = 0.5, 4 steps: the solid
%      stress error is below 1.3e-3, about twice the measured 6.34e-4; with
%      the kinematic row tested with the solid's interface traces it is
%      6.56e-3.
%   6. the snapped circle moving (A = 0.05), crossed grid N = 32, k = 1,
%      dt = 0.025, 8 steps, smooth data: the worst solid stress error over
%      the steps is below 6e-4, about twice the measured 3.17e-4.
%   7. the snapped circle moving (A = 0.1), crossed grid N = 16, k = 1,
%      smooth data, T = 0.2 in 8 and 32 steps (dt = 0.4 h and 0.1 h): the
%      solid stress and displacement errors with 32 steps are below 1.3
%      times those with 8 (measured ratios 0.858, 1.02); with the one-sided
%      differences of u* at every dt ("ustar_bdf") they are 0.908, 1.03.
clc; clear;
warning("off", "MATLAB:nearlySingularMatrix");   % trace extension on sliver patches

n_fail = 0;

%% 1. remap exactness on P_k
flat = {FlatCurveHorizontal(0.45), FlatCurveHorizontal(0.55)};
background = StructuredSimplicialGrid(8, 8, "left");
for k = 1:2
    for snap = [false, true]
        for carry = unique([false, snap])
            [d, switched] = remap_error(k, background, flat, [0; 0.1], snap, carry);
            fprintf("1. flat k = %d, snap %d, carried %d: max difference %.2e; switched on/off %d/%d\n", ...
                k, snap, carry, d, switched);
            n_fail = n_fail + check(d <= 1e-10 && isequal(switched, [16, 0]), ...
                sprintf("flat remap k = %d, snap %d, carried %d", k, snap, carry));
        end
    end
end
circle = {CircularCurve([0.5; 0.5], 0.2), CircularCurve([0.53; 0.5], 0.2)};
[d, switched] = remap_error(2, StructuredSimplicialGrid(16, 16, "crossed"), circle, [0.03; 0], true, true);
fprintf("1. circle k = 2, snapped: max difference %.2e; switched on/off %d/%d\n", d, switched);
n_fail = n_fail + check(d <= 1e-10 && sum(switched) > 0, "circle remap");
circle = {CircularCurve([0.5; 0.5], 0.2), CircularCurve([0.55; 0.5], 0.2)};
[d, switched] = remap_error(2, StructuredSimplicialGrid(16, 16, "crossed"), circle, [0.05; 0], true, true);
fprintf("1. circle k = 2, shift 0.05, snapped: max difference %.2e; switched on/off %d/%d\n", d, switched);
n_fail = n_fail + check(d <= 1e-10 && switched(1) > 0, "circle remap with elements switching on");

%% 2. A = 0 against one fixed carve
k = 1;  dt = 0.05;  n_steps = 4;
uu = "sin(pi*x)*sin(pi*y)*cos(t)";  pp = "cos(pi*x)*(y - 0.5)*(1 + t)";
us = ["0.1*sin(2*x + y)*cos(t)"; "0.1*cos(x - y)*(1 + t^2)"];
% the fixed-carve reference is the toolbox trapezoidal Newmark step with its
% own history (the consistent acceleration differs by 1.2e-6)
out = eulerian_solid_march(k, background, 0, dt, n_steps, uu, pp, us, solid_bossak_rho = 1, ...
    history_acceleration = "translated", history_translate = "library", skeleton_history = "library");
[fs, ss] = fixed_carve_march(k, background, dt, n_steps, uu, pp, us);
d_fix = max(abs([dofs(out.fs) - dofs(fs); dofs(out.ss) - dofs(ss)]));
fprintf("2. A = 0 vs one fixed carve: max dof difference %.2e\n", d_fix);
n_fail = n_fail + check(d_fix <= 1e-12, "A = 0 reproduces the fixed-carve march");

%% 3. exact polynomial solution, both sides re-carved
k = 2;  A = 0.1;  dt = 0.05;  n_steps = 6;
d_t = sprintf("%.17g*sin(2*pi*t)", A);
g = "(1 + t + t^2/2)";
us_poly = ["x*(y - " + d_t + " - 0.45)*" + g; "(x^2 - (y - " + d_t + "))*" + g];
for snap = [true, false]
    [~, out] = evalc("eulerian_solid_march(k, background, A, dt, n_steps, ""x*y*(1 + t)"", " + ...
        """(x + y)*(1 + t)"", us_poly, snap = snap, verbose = true)");   % per-step errors, unprinted
    worst = out.errors.';
    fprintf("3. polynomial, snap %d: worst L2 error per step %s; switched (F on/off, S on/off) %s\n", ...
        snap, mat2str(worst, 3), mat2str(sum(out.history(:, 4:7), 1)));
    n_fail = n_fail + check(max(worst) <= 1e-9, sprintf("polynomial march exact, snap %d", snap));
    n_fail = n_fail + check(all(sum(out.history(:, [5 6]), 1) > 0), ...
        sprintf("polynomial march, snap %d: no row switched on both sides", snap));
end

%% 4. smooth data, bounded and at the level of the fixed carve
k = 1;  dt = 0.05;  n_steps = 9;
[~, out] = evalc("eulerian_solid_march(k, background, 0.1, dt, n_steps, uu, pp, us, verbose = true)");
[~, ref] = evalc("eulerian_solid_march(k, background, 0, dt, n_steps, uu, pp, us, verbose = true)");
worst = max(out.step_errors, [], 1);
ratio = worst ./ max(ref.step_errors, [], 1);
bound = [4, 1.65, 2.5, 1.15];
fprintf("4. smooth, A = 0.1, T = %g, worst over steps: |sF-s| %.3e |uF-u| %.3e |sS-s| %.3e |uS-u| %.3e\n", ...
    n_steps*dt, worst);
fprintf("   / A = 0: %s (bounds %s)\n", mat2str(ratio, 3), mat2str(bound));
fprintf("   max |u_h| (fluid, solid) %.3g %.3g; switched (F on/off, S on/off) %s\n", ...
    max(out.history(:, 2:3), [], 1), mat2str(sum(out.history(:, 4:7), 1)));
n_fail = n_fail + check(all(isfinite([out.eF, out.eS])) && max(out.history(:, 2:3), [], "all") < 10, ...
    "smooth: the march is unbounded");
n_fail = n_fail + check(all(sum(out.history(:, 4:7), 1) > 0), "smooth: no switching on some side");
n_fail = n_fail + check(all(ratio < bound), ...
    "smooth: worst error over the steps of (sF, uF, sS, uS) above its bound times the fixed carve's");

%% 5. circle at rest, time-independent data: the solid stress error is small
[~, out] = evalc("eulerian_solid_march(1, StructuredSimplicialGrid(32, 32, ""crossed""), 0, 0.5, 4, " + ...
    """sin(pi*x)*sin(pi*y)"", ""cos(pi*x)*(y - 0.5)"", [""0.1*sin(2*x + y)""; ""0.1*cos(x - y)""], " + ...
    "geometry = ""circle"")");
fprintf("5. circle at rest, N = 32, T = 2: |sS-s| %.3e (bound 1.3e-3)\n", out.eS(1));
n_fail = n_fail + check(out.eS(1) < 1.3e-3, "circle at rest: solid stress error above its bound");

%% 6. moving circle, smooth data: the worst solid stress error is small
[~, out] = evalc("eulerian_solid_march(1, StructuredSimplicialGrid(32, 32, ""crossed""), 0.05, 0.025, 8, " + ...
    "uu, pp, us, geometry = ""circle"", verbose = true)");
fprintf("6. moving circle, N = 32, A = 0.05: worst |sS-s| %.3e (bound 6e-4)\n", max(out.step_errors(:, 3)));
n_fail = n_fail + check(max(out.step_errors(:, 3)) < 6e-4, "moving circle: solid stress error above its bound");

%% 7. moving circle, small time steps: the solid errors do not grow as dt -> 0
grid = StructuredSimplicialGrid(16, 16, "crossed");
e = zeros(2, 2);
for i = 1:2
    n = 8*4^(i - 1);
    out = eulerian_solid_march(1, grid, 0.1, 0.2/n, n, uu, pp, us, geometry = "circle");
    e(i, :) = out.eS(1:2);
end
fprintf("7. moving circle, N = 16, T = 0.2, 8 -> 32 steps: |sS-s| %.3e -> %.3e, |uS-u| %.3e -> %.3e\n", e);
fprintf("   ratios %s (bound 1.3)\n", mat2str(e(2, :)./e(1, :), 3));
n_fail = n_fail + check(all(e(2, :) <= 1.3*e(1, :)), "small dt: solid error grows as dt decreases");

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_fsi_eulerian_solid_march:failed", "%d check(s) failed", n_fail);
end

function [d, switched] = remap_error(k, background, curves, shift, snap, carry)
% translate_state + remap_solid_state of the projections of a P_k field
% at the three Newmark levels, against the projections of its translate
solid = LinearElasticityProblem(k, young_modulus = 10, poisson_ratio = 0.3, ...
    time_scheme = "newmark", density = 1);
p = {@(x, y) x.^2 - 2*x.*y + y; @(x, y) y.^2 + x.*y - x};
if k == 1, p = {@(x, y) x - 2*y + 1; @(x, y) 3*y + x}; end
q = {@(x, y) 2*x - y; @(x, y) x + 0.5*y - 1};
translated = @(f) cellfun(@(h) @(x, y) h(x - shift(1), y - shift(2)), f, UniformOutput = false);
[old, old_keys, old_faces] = carve(background, curves{1}, snap);
[new, new_keys, new_faces] = carve(background, curves{2}, snap);
state.solid = levels(solid, old, k, p, q);
untranslated = state.solid;
state = translate_state(solid, old, state, shift, k);
extra = {};
if carry, extra = {"untranslated", untranslated, "shift", shift}; end
[s, switched] = remap_solid_state(solid, state.solid, old, new, k, old_keys, new_keys, ...
    old_faces, new_faces, extra{:});
expected = levels(solid, new, k, translated(p), translated(q));
names = string(fieldnames(expected)).';
d = max(arrayfun(@(n) max(abs(s.(n) - expected.(n)), [], "all"), names));
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

function [mesh, keys, face_keys] = carve(background, curve, snap)
% the solid side and its unsnapped element and face keys, as
% eulerian_solid_march carves it
if snap
    [~, mesh] = background.split(curve, max_normal_deviation = 90);
else
    [~, mesh] = background.split(curve);
end
keys = squeeze(mean(reshape(mesh.coordinates(:, mesh.elements), 2, 3, []), 2)).';
face_keys = ((mesh.coordinates(:, mesh.faces(1, :)) + mesh.coordinates(:, mesh.faces(2, :)))/2).';
if snap
    mesh = snap_to_curve(mesh, curve, normal_sign = -1);
end
end

function [fs, ss] = fixed_carve_march(k, background, dt, n_steps, uu, pp, us)
% the snapped carve of Sigma = {y = 0.45} kept for the whole march; backward
% Euler, then BDF2 (theta = 2/3 from (4 v^n - v^n-1)/3) for the fluid,
% Newmark for the solid
data = @(t) translating_flat_problems(k, 0.45, 0, t, dt, uu, pp, us, 1, 1, 1, 10, 0.3);
quad = {"n_quadpoints", k + 2, "n_quadpoints_path", k + 3};
[fp, sp, src, ~, curve] = data(0);
[fm, sm] = background.split(curve, max_normal_deviation = 90);
fm = snap_to_curve(fm, curve);
sm = snap_to_curve(sm, curve, normal_sign = -1);
sp.Data.time_scheme = "newmark";
fp.Data.time_scheme = "euler";
state = FluidStructureCoupling.initial_state(fp, sp, fm, sm, curve, src{:}, "Order", k, quad{:});
v_previous = [];
for n = 1:n_steps
    [fp, sp, src, ~, curve] = data(n*dt);
    sp.Data.time_scheme = "newmark";
    fp.Data.time_scheme = "cranknicolson";
    if n == 1, fp.Data.time_scheme = "euler"; end
    solve_state = state;
    if n > 1
        fp.Data.theta = 2/3;
        solve_state.fluid.velocity = (4*state.fluid.velocity - v_previous)/3;
        solve_state.fluid.acceleration = zeros(size(state.fluid.velocity));
    end
    [fs, ss, new_state] = FluidStructureCoupling.solve(fp, sp, fm, sm, curve, src{:}, ...
        "newmark", "quiet", true, "state", solve_state, "patch_weight", "jacobian", quad{:});
    v_previous = state.fluid.velocity;
    state = new_state;
end
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
