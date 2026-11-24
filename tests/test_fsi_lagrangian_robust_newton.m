% test_fsi_lagrangian_robust_newton.m
% The safeguards of the coupled Newton loop of the hyperelastic Lagrangian
% march (lagrangian_carved_march): line_search (damped Newton, det F > 0
% at every quadrature point and no residual growth) and adaptive_time_step
% (a failed step retried with dt/2). Flat clamped Sigma_0, fluid carved
% from a "left" background grid and snapped, k = 1, N = 8, neo-Hookean
% solid with tau = 4/h (below the march's default 10 mu/h, under which
% the large steps of 2, 3 and 6 converge plainly), data from
% lagrangian_march_problems (t^4 factor).
%   1. both safeguards on, on the march of test_fsi_lagrangian_hyperelastic
%      (A = 0.3, dt = 0.1, T = 0.5): every alpha is 1 and no step is
%      rejected, and the march is bit-identical (isequal) to the plain one;
%   2. A = 0.6 in a single step dt = T = 0.5: plain Newton's first update
%      inverts the solid (HyperelasticProblem:inverted); the damped loop
%      takes alpha < 1, converges, and ends with full steps;
%   3. A = 0.5, dt = T = 0.5: the damped loop finds no admissible step
%      (lagrangian_carved_march:stepFailed); with adaptive_time_step the
%      march halves the step, converges in every accepted step, ends at T,
%      and its error is that of the uniform march at its finest step;
%   4. a forced step failure (newton_max_iterations = 1 on the first
%      attempt at a given t), A = 0.3, T = 0.5:
%      a. at t = 0.1 with dt = 0.1 and no recovery: the march is the
%         uniform dt = 0.05 march, bit for bit (Newmark and Crank-Nicolson
%         carry accelerations, the predictor uses the actual step);
%      b. at t = 0.3 with dt = 0.1, recovery after 2 steps (dt 0.1, 0.1,
%         0.05 x 2, 0.1, ...): the solid against a dt = 0.025 reference
%         at most the dt = 0.1 difference, and the curve (the predictor
%         of the last step, a dt = 0.1 step in both) within 5 % of it;
%      c. the same with a Navier-Stokes fluid (A = 0.03, T = 0.6, the beta
%         extrapolation over a step change);
%      b and c with Sigma_h from u_h (interface_geometry = "displacement");
%      d. at t = 0.3, recovery after 1 step: dt doubles only at t = 0.3,
%         a multiple of the doubled step (dt 0.1, 0.1, 0.05, 0.05, 0.1,
%         0.1), and the snapshots fall on the multiples of dt and carry
%         the step's dt, Newton report, snap count and solved curve;
%   5. a step failing at every dt (newton_max_iterations = 1 on every
%      attempt ending in (0.2, 0.3]) down to min_time_step = dt/4: the
%      error lagrangian_carved_march:timeStepFloor names the reason, the dt
%      reached, the floor and the halvings;
%   6. line_search without adaptive_time_step, A = 0.6, dt = T = 0.5,
%      newton_max_iterations = 2: the second iterate is damped, an error
%      (lagrangian_carved_march:dampedNotConverged);
%   7. the linear march of test_fsi_lagrangian_march (dt = 0.1): no vertex
%      unsnapped, fluid paths below one element height (out.snap); a
%      max_path_ratio below the observed ratio stops it.
clc; clear;
warning("off", "MATLAB:nearlySingularMatrix");   % trace extension on sliver patches
warning("off", "HyperelasticProblem:newtonNotConverged");

n_fail = 0;
k = 1;  N = 8;
background = StructuredSimplicialGrid(N, N, "left");
solid_mesh = apply_tag(StructuredSimplicialGrid(1.5*N, 0.75*N, 0, 1, 0, 0.5), ...
    "Interface", @(P) abs(P(2, :) - 0.5) < 1e-14);
march = @(A, dt, T, varargin) lagrangian_carved_march(k, "flat", background, solid_mesh, A, dt, ...
    round(T/dt), "snap", true, "solid_model", "hyperelastic", "material", "neo_hookean", ...
    "solid_tau", 4*2.6/10, varargin{:});      % tau = 4/h: solid_tau in units of mu = 10/2.6
safeguards = {"line_search", true, "adaptive_time_step", true};

%% 1. alpha = 1 is the plain march, bit for bit
plain = march(0.3, 0.1, 0.5);
guarded = march(0.3, 0.1, 0.5, safeguards{:});
alpha = [guarded.newton.alpha];
same = isequal(fingerprint(plain), fingerprint(guarded));
fprintf("1. safeguards on, A = 0.3, dt = 0.1: alphas %s, rejected %d, isequal %d\n", ...
    sprintf("%g ", unique(alpha)), size(guarded.rejected, 1), same);
n_fail = n_fail + check(all(alpha == 1) && isempty(guarded.rejected), "the safeguards acted");
n_fail = n_fail + check(same, "the safeguarded march differs from the plain one");

%% 2. damping rescues an inverting Newton update
fprintf("2. A = 0.6, dt = 0.5, plain: ");
plain_inverts = throws(@() march(0.6, 0.5, 0.5), "HyperelasticProblem:inverted");
fprintf("throws HyperelasticProblem:inverted %d\n", plain_inverts);
damped = march(0.6, 0.5, 0.5, "line_search", true);
report = damped.newton;
fprintf("   damped (%s): r %s; alpha %s; |u_S - u| %.2e\n", report.stop_reason, ...
    sprintf("%.1e ", report.residual), sprintf("%g ", report.alpha), damped.eS(2));
n_fail = n_fail + check(plain_inverts, "plain Newton did not invert");
n_fail = n_fail + check(report.converged && report.stop_reason == "tolerance" ...
    && min(report.alpha) < 1 && report.alpha(end) == 1, "damped Newton did not converge with damping");

%% 3. halving the step where damping alone fails
fprintf("3. A = 0.5, dt = 0.5, damped: ");
damped_fails = throws(@() march(0.5, 0.5, 0.5, "line_search", true), "lagrangian_carved_march:stepFailed");
fprintf("throws stepFailed %d\n", damped_fails);
adaptive = march(0.5, 0.5, 0.5, safeguards{:});
dts = adaptive.time_steps.';
uniform = march(0.5, min(dts), 0.5);
fprintf("   adaptive: rejected %d, dt %s, t_end %.17g; |u_S - u| %.3e (uniform dt = %g: %.3e)\n", ...
    size(adaptive.rejected, 1), sprintf("%g ", dts), adaptive.history(end, 1), ...
    adaptive.eS(2), min(dts), uniform.eS(2));
n_fail = n_fail + check(damped_fails, "damped Newton alone did not fail");
n_fail = n_fail + check(~isempty(adaptive.rejected) && all([adaptive.newton.converged]) ...
    && adaptive.history(end, 1) == 0.5, "the adaptive march did not recover");
n_fail = n_fail + check(adaptive.eS(2) < 1.5*uniform.eS(2), "adaptive error above the uniform finest-step level");

%% 4a. a forced halving at the first step, no recovery = the uniform dt/2 march
forced = march(0.3, 0.1, 0.5, "problems", forced_failure(0.3, 0.1, 0.1, "stokes"), ...
    "adaptive_time_step", true, "time_step_recovery", 1000);
half = march(0.3, 0.05, 0.5);
same = isequal(fingerprint(forced), fingerprint(half));
fprintf("4a. forced at t = 0.1, no recovery: dt %s; isequal to uniform dt = 0.05: %d\n", ...
    sprintf("%g ", forced.time_steps), same);
n_fail = n_fail + check(size(forced.rejected, 1) == 1 && all(forced.time_steps == 0.05), ...
    "the forced failure was not retried at dt/2");
n_fail = n_fail + check(same, "the halved march differs from the uniform dt/2 march");

%% 4b, 4c. a forced halving mid-march, then recovery
for c = {{0.3, 0.5, "stokes"}, {0.03, 0.6, "navier_stokes"}}
    [A, T, fluid] = c{1}{:};
    run = @(dt, varargin) march(A, dt, T, "fluid_model", fluid, "interface_geometry", "displacement", varargin{:});
    reference = run(0.025);
    distance = @(r) [max(abs(r.ss.Unknowns.field(2).dof - reference.ss.Unknowns.field(2).dof), [], "all"), ...
        max(vecnorm(r.curve.get_2d_points(linspace(0, 1, 201)) ...
        - reference.curve.get_2d_points(linspace(0, 1, 201))))];
    coarse = distance(run(0.1));
    fine = distance(run(0.05));
    forced = run(0.1, "problems", forced_failure(A, 0.1, 0.3, fluid), ...
        "adaptive_time_step", true, "time_step_recovery", 2);
    shifted = distance(forced);
    fprintf("4%s. %s, forced at t = 0.3: dt %s\n", char('b' + (fluid ~= "stokes")), fluid, ...
        sprintf("%g ", forced.time_steps));
    fprintf("   |u_S - u_ref|, |Sigma - Sigma_ref|: dt 0.1 %s, forced %s, dt 0.05 %s\n", ...
        sprintf("%.2e ", coarse), sprintf("%.2e ", shifted), sprintf("%.2e ", fine));
    n_fail = n_fail + check(isequal(forced.time_steps.', [0.1 0.1 0.05 0.05 0.1 0.1*ones(1, round(T/0.1) - 4)]) ...
        && abs(forced.history(end, 1) - T) < 1e-14, fluid + ": unexpected step sequence");
    n_fail = n_fail + check(shifted(1) <= coarse(1) && shifted(2) <= 1.05*coarse(2), ...
        fluid + ": a step change raised the error");
end

%% 4d. recovery waits for the grid: dt doubles only at a multiple of the doubled step
snapshot_file = string(tempname) + ".mat";
forced = march(0.3, 0.1, 0.5, "problems", forced_failure(0.3, 0.1, 0.3, "stokes"), ...
    "adaptive_time_step", true, "time_step_recovery", 1, "snapshot_file", snapshot_file);
saved = load(snapshot_file);
delete(snapshot_file);
entries = cellfun(@(name) saved.(name), sort(fieldnames(saved)));
saved_t = [entries.t];
fprintf("4d. forced at t = 0.3, recovery after 1 step: dt %s; snapshots at t = %s\n", ...
    sprintf("%g ", forced.time_steps), sprintf("%g ", saved_t));
n_fail = n_fail + check(isequal(forced.time_steps.', [0.1 0.1 0.05 0.05 0.1 0.1]), ...
    "recovery left the grid of dt");
n_fail = n_fail + check(numel(saved_t) == 5 && max(abs(saved_t - (1:5)*0.1)) < 1e-14, ...
    "snapshots off the multiples of dt");
fprintf("    snapshot dt %s, Newton iterations %s, stop reasons %s\n", sprintf("%g ", [entries.dt]), ...
    sprintf("%d ", [entries.newton_iterations]), strjoin(string({entries.stop_reason})));
n_fail = n_fail + check(isequal([entries.dt], [0.1 0.1 0.05 0.1 0.1]) ...
    && isequal([entries.newton_iterations], [forced.newton([1 2 4 5 6]).iterations]) ...
    && all(string({entries.stop_reason}) == [forced.newton([1 2 4 5 6]).stop_reason]) ...
    && isequal([entries.unsnapped], forced.snap.unsnapped([1 2 4 5 6]).') ...
    && isequal(size(entries(end).solved_sigma_points), size(entries(end).sigma_points)) ...
    && max(vecnorm(entries(end).solved_sigma_points - forced.curve.get_2d_points(linspace(0, 1, 401)))) > 0, ...
    "snapshot diagnostics");

%% 5. a step failing down to min_time_step is an error naming the reason
caught = error_of(@() march(0.3, 0.1, 0.5, "problems", forced_failure(0.3, 0.1, 0.3, "stokes", true), ...
    "adaptive_time_step", true, "min_time_step", 0.025));
fprintf("5. failing at every dt: %s\n   %s\n", caught.identifier, caught.message);
n_fail = n_fail + check(caught.identifier == "lagrangian_carved_march:timeStepFloor" ...
    && contains(caught.message, "(max_iterations) at dt = 0.025 = dt/2^2") ...
    && contains(caught.message, "floor min_time_step = 0.025") && contains(caught.message, "2 halving(s)"), ...
    "the floor error does not name reason, dt, floor and halvings");

%% 6. damped and out of iterations, no adaptive step: an error, not a damped state
caught = error_of(@() march(0.6, 0.5, 0.5, "line_search", true, ...
    "problems", forced_failure(0.6, 0.5, 0.5, "stokes", false, 2)));
fprintf("6. A = 0.6, damped, newton_max_iterations = 2: %s\n   %s\n", caught.identifier, caught.message);
n_fail = n_fail + check(caught.identifier == "lagrangian_carved_march:dampedNotConverged", ...
    "a damped step out of iterations was accepted");

%% 7. snap diagnostics and max_path_ratio
linear = lagrangian_carved_march(k, "flat", background, solid_mesh, 0.03, 0.1, 5, "snap", true);
ratios = linear.snap.path_ratio.';
fprintf("7. linear march: unsnapped %s; path/h_perp %s\n", sprintf("%d ", linear.snap.unsnapped), ...
    sprintf("%.2e ", ratios));
n_fail = n_fail + check(all(linear.snap.unsnapped == 0) && numel(ratios) == 5 && all(ratios >= 0 & ratios < 1) ...
    && max(ratios) > 0, ...
    "unexpected snap diagnostics");
caught = error_of(@() lagrangian_carved_march(k, "flat", background, solid_mesh, 0.03, 0.1, 5, "snap", true, ...
    "max_path_ratio", max(ratios)/2));
fprintf("   max_path_ratio = %.2e: %s\n", max(ratios)/2, caught.message);
n_fail = n_fail + check(caught.identifier == "lagrangian_carved_march:pathRatio", "max_path_ratio did not stop the march");

fprintf("\n%d check(s) failed\n", n_fail);
if n_fail > 0
    error("test_fsi_lagrangian_robust_newton:failed", "%d check(s) failed", n_fail);
end

function err = error_of(f)
% the exception f() raises (an empty-identifier one when it raises none)
err = MException("", "no error");
try
    f();
catch err
end
end

function f = fingerprint(r)
% the numeric output of a march
f = {r.eF, r.eS, r.paths, r.history, {r.newton.residual}, ...
    arrayfun(@(x) x.dof, r.fs.Unknowns.field, "UniformOutput", false), ...
    arrayfun(@(x) x.dof, r.ss.Unknowns.field, "UniformOutput", false), r.ss.Unknowns.hybrid(1).dof};
end

function problems = forced_failure(A, dt, t_fail, fluid, every, iterations)
% the march data, with newton_max_iterations = iterations on the first
% call at t_fail, or (every = true) on every call at t in
% (t_fail - dt, t_fail]
arguments
    A, dt, t_fail, fluid
    every (1,1) logical = false
    iterations (1,1) double = 1
end
called = containers.Map("KeyType", "double", "ValueType", "logical");   % a handle: shared by the calls
problems = @(t) forced_data(t, A, dt, t_fail, fluid, called, every, iterations);
end

function [fp, sp, pd, reference, sources] = forced_data(t, A, dt, t_fail, fluid, called, every, iterations)
[fp, sp, pd, reference, sources] = lagrangian_march_problems(1, "flat", t, dt, A, ...
    "solid_model", "hyperelastic", "material", "neo_hookean", ...
    "convection", fluid == "navier_stokes");
if every && t > t_fail - dt + 1e-12 && t < t_fail + 1e-12
    sp.Data.newton_max_iterations = iterations;
elseif abs(t - t_fail) < 1e-12 && ~isKey(called, t_fail)
    called(t_fail) = true; %#ok<NASGU> (a handle)
    sp.Data.newton_max_iterations = iterations;
end
end
