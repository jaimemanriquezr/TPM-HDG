% test_biot_dynamic.m
% Fully dynamic Biot (BiotProblem with density = rho_P > 0): the inertia
% rho_P d_tt eta is the backward difference applied twice,
%   v^{n+1} = (eta^{n+1} - eta*)/cdt,  a^{n+1} = (v^{n+1} - v*)/cdt,
% with (cdt, eta*, v*) from biot_step_history.
%
%   (a) biot_step_history combines the volume displacement and velocity;
%       initial_state and get_kinematic_state carry them;
%   (b) temporal order (V1): data of degree <= k in space, smooth in time,
%       rho_P = 1, T = 1, dt = 1/4 ... 1/32, k = 2, 8 x 8: EOC of eta, v
%       and p between the two finest steps >= 0.9 (euler), >= 1.8 (bdf2
%       from an exact two-level start), >= 1.2 (bdf2 with its Euler start);
%   (c) spatial order (V2): one Euler step from the exact state with data
%       linear in t (the backward differences are then exact), h = 1/4,
%       1/8, 1/16, k = 1, 2, 3: rate >= k + 1 - 0.3 (eta, p),
%       >= k + 1 - 0.5 (sigma, q);
%   (d) energy: zero data, nonzero initial displacement, all-Dirichlet box,
%       euler, 20 steps: 1/2 rho_P ||v||^2 + 1/2 (A sigma_E, sigma_E)
%       + 1/2 s0 ||p||^2 is nonincreasing.
% Unit square, all-Dirichlet, alpha = 1, kappa = I.
clc; clear;
n_fail = 0;

%% (a) bookkeeping
B = struct("displacement", 3*ones(2), "displacement_previous", ones(2), ...
    "velocity", 2*ones(2), "velocity_previous", 5*ones(2), ...
    "displacement_skeleton_previous", 1);
[cdt, h] = biot_step_history("bdf2", 0.3, B);
ok = abs(cdt - 0.2) < 1e-15 ...
    && max(abs(h.displacement - (11/3)*ones(2)), [], "all") < 1e-14 ...
    && max(abs(h.velocity - 1*ones(2)), [], "all") < 1e-14;
n_fail = n_fail + check(ok, "bdf2: displacement* and velocity* = (4 w^n - w^{n-1})/3");
[~, h] = biot_step_history("euler", 0.3, B);
ok = isequal(h.displacement, 3*ones(2)) && isequal(h.velocity, 2*ones(2));
n_fail = n_fail + check(ok, "euler: displacement* = displacement, velocity* = velocity");
[~, h] = biot_step_history("euler", 0.3);
ok = isequal(h.displacement, 0) && isequal(h.velocity, 0);
n_fail = n_fail + check(ok, "empty state: zero displacement and velocity history");

k = 2; N = 8; dt = 0.25;
E = 1; nu = 0.3; alpha = 1; s0 = 1; kappa = eye(2); rho = 1;
params = {"young_modulus", E, "poisson_ratio", nu, "biot_alpha", alpha, ...
    "storage", s0, "diffusion", kappa};
mesh = StructuredSimplicialGrid(N, N);       % every boundary face "Dirichlet"
moving = BiotProblem.get_manufactured_solution( ...
    ["exp(t)*(x^2 + x*y)"; "exp(t)*(x*y - y^2 + x)"], ...
    "cos(2*t)*(x^2 - x*y + y^2 + 1)", params{:}, density = rho);
for scheme = ["euler", "bdf2"]
    problem = BiotProblem(k, params{:}, density = rho, time_step = dt, time_scheme = scheme);
    problem = set_initial(problem, moving.at(0), 0);
    state0 = initial_state(problem, mesh, k);
    b0 = state0.biot;
    d = moving.at(0);
    ok = max(abs(b0.displacement - L2proj2D(mesh, d.u_P, k)), [], "all") < 1e-14 ...
        && max(abs(b0.velocity - L2proj2D(mesh, d.v_P, k)), [], "all") < 1e-14 ...
        && isempty(b0.displacement_previous) && isempty(b0.velocity_previous);
    n_fail = n_fail + check(ok, scheme + ": initial_state displacement and velocity");
    state = state0;
    for n = 1:2
        problem = set_sources(problem, moving.at(state.time + dt));
        problem.Data.state = state;
        solution = solve(problem, mesh, k);
        [c, h] = biot_step_history(scheme, dt, state.biot);
        previous = state.biot;
        state = get_kinematic_state(problem, solution, state);
        eta1 = solution.Unknowns.field(2).dof;
        b = state.biot;
        ok = isequal(b.displacement, eta1) ...
            && max(abs(b.velocity - (eta1 - h.displacement)/c), [], "all") < 1e-12 ...
            && isequal(b.displacement_previous, previous.displacement) ...
            && isequal(b.velocity_previous, previous.velocity);
        n_fail = n_fail + check(ok, sprintf("%s, step %d: displacement, velocity = (eta - eta*)/cdt", ...
            scheme, n));
    end
end

%% (b) temporal order (V1)
% "bdf2" starts with its Euler step (the toolbox's start); "bdf2, exact
% start" starts from the projected and lifted states at t = 0 and t = dt
% and isolates the order of the BDF2 march. The Euler start costs order in
% this range of dt (v: EOC about 1.3, p: about 1.7 between 1/16 and 1/32),
% so its bound is 1.2, the plan's 1.8 is asserted on the exact start.
T = 1;
dts = 1./[4, 8, 16, 32];
thresholds = struct("euler", 0.9, "bdf2", 1.2, "bdf2_exact_start", 1.8);
labels = ["eta", "v", "p"];
history_names = ["displacement_skeleton", "stress", "rotation", "pressure", ...
    "displacement", "velocity"];
for variant = ["euler", "bdf2", "bdf2_exact_start"]
    scheme = extractBefore(variant + "_", "_");
    errors = zeros(numel(dts), 3);
    for i = 1:numel(dts)
        dt = dts(i);
        problem = BiotProblem(k, params{:}, density = rho, time_step = dt, time_scheme = scheme);
        problem = set_initial(problem, moving.at(0), 0);
        state = initial_state(problem, mesh, k);
        if variant == "bdf2_exact_start"
            problem = set_initial(problem, moving.at(dt), dt);
            previous = state;
            state = initial_state(problem, mesh, k);
            state.step = 1;
            for name = history_names
                state.biot.(name + "_previous") = previous.biot.(name);
            end
        end
        for n = state.step + 1:round(T/dt)
            problem = set_sources(problem, moving.at(state.time + dt));
            problem.Data.state = state;
            solution = solve(problem, mesh, k);
            state = get_kinematic_state(problem, solution, state);
        end
        d = moving.at(state.time);
        velocity = solution.Unknowns.field(2);
        e_eta = compute_l2_error(velocity, mesh, Datum(d.u_P));
        velocity.dof = state.biot.velocity;
        errors(i, :) = [e_eta, compute_l2_error(velocity, mesh, Datum(d.v_P)), ...
            compute_l2_error(solution.Unknowns.field(5), mesh, Datum(d.p_P))];
    end
    eoc = log2(errors(1:end-1, :)./errors(2:end, :));
    fprintf("\n(b) %s: errors at T = %g, rho_P = %g\n%8s", strrep(variant, "_", " "), T, rho, "dt");
    fprintf(" %12s %6s", [labels; repmat("EOC", 1, 3)]);
    fprintf("\n");
    for i = 1:numel(dts)
        fprintf("%8.5f", dts(i));
        for j = 1:3
            if i == 1
                fprintf(" %12.4e %6s", errors(i, j), "-");
            else
                fprintf(" %12.4e %6.2f", errors(i, j), eoc(i-1, j));
            end
        end
        fprintf("\n");
    end
    threshold = thresholds.(variant);
    n_fail = n_fail + check(all(eoc(end, :) >= threshold), ...
        sprintf("(b) %s: EOC in time of eta, v, p (two finest steps) >= %.1f", ...
        strrep(variant, "_", " "), threshold));
end

%% (c) spatial order (V2)
t0 = 0.5; dt = 0.25;
names = ["sigma", "eta", "q", "p"];
linear = BiotProblem.get_manufactured_solution( ...
    ["(1 + t)*sin(pi*x)*cos(pi*y)"; "(2 - t)*x*exp(y)"], ...
    "(1 + 2*t)*cos(pi*x)*sin(pi*y) + t*x*y", params{:}, density = rho);
levels = [4 8 16];
for k = 1:3
    err = zeros(numel(levels), 4);
    for i = 1:numel(levels)
        mesh = StructuredSimplicialGrid(levels(i), levels(i));
        problem = BiotProblem(k, params{:}, density = rho, time_step = dt);
        problem = set_initial(problem, linear.at(t0), t0);
        problem.Data.state = initial_state(problem, mesh, k);
        d = linear.at(t0 + dt);
        problem = set_sources(problem, d);
        sol = solve(problem, mesh, k);
        err(i, :) = [compute_l2_error(sol.Unknowns.field(1), mesh, Datum(d.sigma_P)), ...
            compute_l2_error(sol.Unknowns.field(2), mesh, Datum(d.u_P)), ...
            compute_l2_error(sol.Unknowns.field(4), mesh, Datum(d.q_P)), ...
            compute_l2_error(sol.Unknowns.field(5), mesh, Datum(d.p_P))];
    end
    eoc = [nan(1, 4); log2(err(1:end-1, :)./err(2:end, :))];
    fprintf("\n(c) k = %d, one Euler step from the exact state at t = %g\n  %4s  %s\n", k, t0, "n", ...
        strjoin(compose("%10s  %5s", names.', repmat("eoc", 4, 1)).', "  "));
    for i = 1:numel(levels)
        fprintf("  %4d  %s\n", levels(i), ...
            strjoin(compose("%10.3e  %5.2f", err(i, :).', eoc(i, :).').', "  "));
    end
    rate = log(err(1, :)./err(end, :))/log(levels(end)/levels(1));
    fprintf("  rate coarsest-finest: %s\n", strjoin(names + " " + compose("%.2f", rate), ", "));
    ok = all(rate([2 4]) >= k + 1 - 0.3) && all(rate([1 3]) >= k + 1 - 0.5);
    n_fail = n_fail + check(ok, "(c) rates k = " + k);
end

%% (d) energy
k = 2; N = 8; dt = 0.05; n_steps = 20;
mesh = StructuredSimplicialGrid(N, N);
problem = BiotProblem(k, params{:}, density = rho, time_step = dt);
problem.Data.initial_displacement = {@(x, y) sin(pi*x).*sin(pi*y); ...
    @(x, y) x.*(1 - x).*y.*(1 - y)};
state = initial_state(problem, mesh, k);
energy = zeros(1, n_steps + 1);
energy(1) = biot_energy(state.biot, mesh, k, E, nu, alpha, s0, rho);
for n = 1:n_steps
    problem.Data.state = state;
    solution = solve(problem, mesh, k);
    state = get_kinematic_state(problem, solution, state);
    energy(n + 1) = biot_energy(state.biot, mesh, k, E, nu, alpha, s0, rho);
end
increase = max(diff(energy))/energy(1);
fprintf("\n(d) energy: E^0 = %.6e, E^%d = %.6e, max relative increase per step = %.3e\n", ...
    energy(1), n_steps, energy(end), increase);
n_fail = n_fail + check(increase <= 1e-12, "(d) energy nonincreasing over 20 Euler steps");

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

function problem = set_initial(problem, d, t)
% the initial data from the manufactured fields at time t
problem.Data.time = t;
problem.Data.initial_displacement = d.u_P;
problem.Data.initial_velocity = d.v_P;
problem.Data.initial_pressure = d.p_P;
end

function problem = set_sources(problem, d)
% sources and Dirichlet data from the manufactured fields at t_{n+1}
problem.Data.source_term = d.f_P;
problem.Data.source_mass = d.g_P;
problem.Data.datum_dirichlet = [d.u_P; {d.p_P}];
end

function e = biot_energy(b, mesh, k, E, nu, alpha, s0, rho)
% 1/2 rho ||v||^2 + 1/2 (A sigma_E, sigma_E) + 1/2 s0 ||p||^2, with
% sigma_E = sigma + alpha p I and (A ., .) the [1, 1] block of
% LinearElasticityMatrices
M = pages(mass_matrix(mesh, DubinerBasis(k), DubinerBasis(k)));
n_el = size(mesh.elements, 2);
quad = @(A, w) sum(pagemtimes(reshape(w, 1, [], 1, 1, n_el), ...
    pagemtimes(A, reshape(w, [], 1, 1, 1, n_el))), "all");
lhs = LinearElasticityMatrices(k, mesh, E, nu);
A = [];
for i = 1:numel(lhs.Arrays)
    if isequal(lhs.Arrays{i}.id, [1 1]), A = double(lhs.Arrays{i}); end
end
sigma_e = BiotProblem.effective_stress(b.stress, b.pressure, alpha, k);
e = rho/2*quad(pagekron(eye(2), M), b.velocity) + quad(A, sigma_e)/2 ...
    + s0/2*quad(M, b.pressure);
end
