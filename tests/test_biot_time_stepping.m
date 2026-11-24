% test_biot_time_stepping.m
% Time stepping of the quasi-static Biot problem (BiotProblem): the
% history combinations of biot_step_history, the initial state
% (initial_state: projections and lifted total stress) and its advance
% (get_kinematic_state), under backward Euler and BDF2.
%
%   (a) the tiny cases of biot_step_history;
%   (b) steady data: from initial_state, the exact discrete state is kept
%       over 5 steps (both schemes), and get_kinematic_state copies the
%       solved fields and moves the current ones to "_previous";
%   (c) temporal order: data polynomial of degree <= k in space (no
%       spatial error in the exact state), smooth in time, T = 1,
%       Δt = 1/4 ... 1/32: EOC of η and p at T between the two finest
%       steps >= 0.9 (euler), >= 1.8 (bdf2); the rate over all steps is
%       printed (Euler is pre-asymptotic at Δt = 1/4).
% Unit square, k = 2, 8 x 8, all-Dirichlet, α = 1, s0 = 1, κ = I.
clc; clear;
n_fail = 0;

%% (a) biot_step_history, tiny cases
[cdt, h] = biot_step_history("bdf2", 0.3);
ok = cdt == 0.3 && all(arrayfun(@(f) isequal(h.(f), 0), ...
    ["displacement_skeleton", "stress", "rotation", "pressure"]));
n_fail = n_fail + check(ok, "bdf2, empty state: cdt = dt, zero history");

B = struct("stress", 3*ones(2), "stress_previous", ones(2), ...
    "displacement_skeleton_previous", 1);
[cdt, h] = biot_step_history("bdf2", 0.3, B);
ok = abs(cdt - 0.2) < 1e-15 && max(abs(h.stress - (11/3)*ones(2)), [], "all") < 1e-14 ...
    && isequal(h.pressure, 0) && isequal(h.rotation, 0) ...
    && abs(h.displacement_skeleton - (-1/3)) < 1e-15;
n_fail = n_fail + check(ok, "bdf2: cdt = 2 dt/3, stress* = (4*3 - 1)/3");

[cdt, h] = biot_step_history("euler", 0.3, B);
ok = cdt == 0.3 && isequal(h.stress, 3*ones(2)) && isequal(h.pressure, 0) ...
    && isequal(h.displacement_skeleton, 0);
n_fail = n_fail + check(ok, "euler: cdt = dt, stress* = stress");

B.displacement_skeleton_previous = [];
[cdt, h] = biot_step_history("bdf2", 0.3, B);
ok = cdt == 0.3 && isequal(h.stress, 3*ones(2));
n_fail = n_fail + check(ok, "bdf2 without level n-1: the Euler step");

%% common setup
k = 2; N = 8;
E = 1; nu = 0.3; alpha = 1; s0 = 1; kappa = eye(2);
mesh = StructuredSimplicialGrid(N, N);       % every boundary face "Dirichlet"
d1 = k + 1;
params = {"young_modulus", E, "poisson_ratio", nu, "biot_alpha", alpha, ...
    "storage", s0, "diffusion", kappa};

%% (b) steady data
steady = biot_data(["x.^2 + x.*y - y"; "x.*y - y.^2 + 2*x"], "x.^2 - x.*y + y.^2 + 1", ...
    E, nu, alpha, s0, kappa);
dt = 0.2;
for scheme = ["euler", "bdf2"]
    problem = BiotProblem(k, params{:}, time_step = dt, time_scheme = scheme);
    problem = set_data(problem, steady, 0, true);
    state0 = initial_state(problem, mesh, k);
    state = state0;
    err = 0;
    for n = 1:5
        problem = set_data(problem, steady, state.time + dt, false);
        problem.Data.state = state;
        solution = solve(problem, mesh, k);
        previous = state;
        state = get_kinematic_state(problem, solution, state);
        err = max([err, state_distance(state.biot, state0.biot), ...
            field_errors(solution, mesh, steady, state.time)]);
        if n == 1
            % the advance copies the solved fields and keeps level n-1
            u = solution.Unknowns;
            b = state.biot;
            ok = isequal(b.stress, u.field(1).dof) && isequal(b.rotation, u.field(3).dof) ...
                && isequal(b.pressure, u.field(5).dof) ...
                && isequal(b.displacement_skeleton, u.hybrid(1).dof(1:2*d1, :)) ...
                && isequal(b.pressure_skeleton, u.hybrid(1).dof(2*d1+1:end, :)) ...
                && isequal(b.stress_previous, previous.biot.stress) ...
                && isequal(b.rotation_previous, previous.biot.rotation) ...
                && isequal(b.pressure_previous, previous.biot.pressure) ...
                && isequal(b.displacement_skeleton_previous, previous.biot.displacement_skeleton) ...
                && state.step == 1 && abs(state.time - dt) < 1e-15;
            n_fail = n_fail + check(ok, scheme + ": get_kinematic_state field copies");
        end
    end
    fprintf("%s, steady data, 5 steps: max dof / L2 deviation = %.3e\n", scheme, err);
    n_fail = n_fail + check(err < 1e-9, scheme + ": steady state kept to 1e-9");
end

% an empty input state is the zero history at Data.time, step 0
problem = BiotProblem(k, params{:}, time_step = dt);
problem = set_data(problem, steady, dt, false);
solution = solve(problem, mesh, k);
state = get_kinematic_state(problem, solution);
ok = state.step == 1 && abs(state.time - dt) < 1e-15 ...
    && isequal(state.biot.stress_previous, zeros(size(state.biot.stress))) ...
    && isequal(state.biot.displacement_skeleton_previous, ...
               zeros(size(state.biot.displacement_skeleton)));
n_fail = n_fail + check(ok, "empty state: step from rest, zero level n-1");

%% (c) temporal order
moving = biot_data(["exp(t).*(x.^2 + x.*y)"; "exp(t).*(x.*y - y.^2 + x)"], ...
    "cos(2*t).*(x.^2 - x.*y + y.^2 + 1)", E, nu, alpha, s0, kappa);
T = 1;
dts = 1./[4, 8, 16, 32];
thresholds = struct("euler", 0.9, "bdf2", 1.8);
for scheme = ["euler", "bdf2"]
    errors = zeros(numel(dts), 2);
    for i = 1:numel(dts)
        dt = dts(i);
        problem = BiotProblem(k, params{:}, time_step = dt, time_scheme = scheme);
        problem = set_data(problem, moving, 0, true);
        state = initial_state(problem, mesh, k);
        for n = 1:round(T/dt)
            problem = set_data(problem, moving, state.time + dt, false);
            problem.Data.state = state;
            solution = solve(problem, mesh, k);
            state = get_kinematic_state(problem, solution, state);
        end
        errors(i, :) = field_errors(solution, mesh, moving, state.time);
    end
    eoc = log2(errors(1:end-1, :)./errors(2:end, :));
    eoc_all = log2(errors(1, :)./errors(end, :))/(numel(dts) - 1);
    fprintf("\n%s: errors at T = %g\n%8s %12s %6s %12s %6s\n", scheme, T, "dt", "eta", "EOC", "p", "EOC");
    for i = 1:numel(dts)
        if i == 1
            fprintf("%8.5f %12.4e %6s %12.4e %6s\n", dts(i), errors(i, 1), "-", errors(i, 2), "-");
        else
            fprintf("%8.5f %12.4e %6.2f %12.4e %6.2f\n", dts(i), errors(i, 1), eoc(i-1, 1), ...
                errors(i, 2), eoc(i-1, 2));
        end
    end
    fprintf("%8s %12s %6.2f %12s %6.2f   (coarsest to finest)\n", "all", "", eoc_all(1), "", eoc_all(2));
    threshold = thresholds.(scheme);
    n_fail = n_fail + check(all(eoc(end, :) >= threshold), ...
        sprintf("%s: EOC in time of eta and p (two finest steps) >= %.1f", scheme, threshold));
end

fprintf("\n%d check(s) failed.\n", n_fail);

%% local functions
function n = check(ok, label)
if ok
    fprintf("  ok    %s\n", label);
    n = 0;
else
    fprintf("  FAIL  %s\n", label);
    n = 1;
end
end

function d = biot_data(displacement, pressure, E, nu, alpha, s0, kappa)
% exact Biot fields and data as handles of (x, y, t): u {2,1}, p, f {2,1}, g
syms x y t real
[lambda, mu] = LinearElasticityProblem.get_lame_parameters(E, nu);
eta = str2sym(displacement);
p = str2sym(pressure);
grad_eta = jacobian(eta, [x, y]);
div_eta = trace(grad_eta);
sigma = mu*(grad_eta + grad_eta.') + lambda*div_eta*eye(2) - alpha*p*eye(2);
f = -[divergence(sigma(1, :), [x, y]); divergence(sigma(2, :), [x, y])];
q = -kappa*gradient(p, [x, y]);
g = s0*diff(p, t) + alpha*diff(div_eta, t) + divergence(q, [x, y]);
to_handle = @(e) broadcast(matlabFunction(e, Vars = [x, y, t]));
d.u = {to_handle(eta(1)); to_handle(eta(2))};
d.p = to_handle(p);
d.f = {to_handle(f(1)); to_handle(f(2))};
d.g = to_handle(g);
end

function h = broadcast(fun)
h = @(x, y, t) fun(x, y, t) + 0*x;
end

function problem = set_data(problem, d, t, initial)
% the data at time t: sources and Dirichlet data, or (initial) the initial data
at = @(h) @(x, y) h(x, y, t);
if initial
    problem.Data.time = t;
    problem.Data.initial_displacement = {at(d.u{1}); at(d.u{2})};
    problem.Data.initial_pressure = at(d.p);
else
    problem.Data.source_term = {at(d.f{1}); at(d.f{2})};
    problem.Data.source_mass = at(d.g);
    problem.Data.datum_dirichlet = {at(d.u{1}); at(d.u{2}); at(d.p)};
end
end

function e = field_errors(solution, mesh, d, t)
% L2 errors of eta and p at time t
at = @(h) @(x, y) h(x, y, t);
e = [compute_l2_error(solution.Unknowns.field(2), mesh, Datum({at(d.u{1}); at(d.u{2})})), ...
     compute_l2_error(solution.Unknowns.field(5), mesh, Datum(at(d.p)))];
end

function d = state_distance(a, b)
% max dof difference of the current fields of two state.biot structs
names = ["displacement_skeleton", "stress", "rotation", "pressure", "pressure_skeleton"];
d = 0;
for name = names
    d = max(d, max(abs(a.(name) - b.(name)), [], "all"));
end
end
