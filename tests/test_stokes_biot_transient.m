% test_stokes_biot_transient.m
% Transient Stokes-Biot coupling (StokesBiotCoupling.solve with fluid
% density rho_F > 0 and Biot density rho_P > 0): the fluid inertia
% rho_F (u - u*)/cdt and the Biot inertia take one scheme (Euler or BDF2)
% and one step for both media, read from the Biot Data; (cdt, u*) come
% from fluid_step_history, (cdt, eta*, v*) from biot_step_history.
% Unit square, Sigma: y = 1/2, fluid above, Biot below, fitted_split
% meshes of sizes (0.8 h, h), Dirichlet data on the whole outer boundary,
% pinned fluid pressure, s0 = 1, friction_bjs = 1, alpha = 1, kappa = I.
%
%   (a) bookkeeping: fluid_step_history (Euler, BDF2, sentinel, zeros,
%       the extrapolated convective velocity); StokesBiotCoupling
%       .initial_state and .advance_state (time and step advanced once,
%       u_star kept); the errors StokesBiotCoupling:timeStep and
%       :historyLevels; V0: with rho_F = 0 the step does not read
%       state.fluid and reproduces the errors recorded on the code before
%       the transient branch (to 1e-10 relative: the recorded numbers are
%       reproduced to about 4e-12 on another machine and MATLAB release);
%   (b) patch test, flat Sigma, no gap, rho_F = rho_P = 1, k = 1, 2, 3,
%       h = 1/4, 1/6: fields polynomial of degree <= k in space (stream
%       function k + 1), every field <= 1e-10 relative:
%         euler, linear in t, 3 steps from the initial state at t = 1/2
%           (exercises the fluid history u^n != 0 and the Biot history
%           eta*, v* with a zero discrete acceleration);
%         bdf2 with its Euler start, linear in t, 3 steps (the first step
%           Euler on both media, then BDF2);
%         bdf2 from an exact two-level start, quadratic in t, 2 steps
%           (BDF2 is exact for quadratics, so d_tt eta = const != 0: the
%           Biot inertia and the fluid inertia are exercised in full);
%   (c) time order (V6 without convection): stream function of degree 3,
%       eta, p_F, p_P of degree 2 in space, smooth in t, gap h^2, k = 2,
%       h = 1/4, T = 1, dt = 1/4 ... 1/128; EOC of u_F, sigma_F, eta, v,
%       sigma_P, q, p between the two finest steps: >= 0.9 (euler),
%       >= 1.8 (bdf2, exact two-level start), reported (bdf2 with the
%       toolbox's Euler start). The range 1/4 ... 1/32 is pre-asymptotic
%       for eta and sigma_P with these data (EOC 1/16 -> 1/32: 0.86, 0.87
%       euler; 1.71 eta bdf2 exact start, also without a gap), hence the
%       two finer steps;
%   (d) V5 without convection: 3 bdf2 steps from an exact two-level start
%       with exact time-derivative data (smooth in space, quadratic in t,
%       for which BDF2 is exact), gap h^2, h = 1/8, k = 2, dt = 1/20: the
%       error of every field after 3 steps is within 2x the error after
%       the first step (the toolbox's Euler start reported). The start
%       state is the projection of the exact fields, not the discrete
%       solution: their O(h^(k+1)) difference enters v as about 1/cdt
%       times it, so at much smaller dt the v error grows over the first
%       steps before it saturates (a start transient, not amplification);
%   (e) sanity space rate (reported): one Euler step from the exact state,
%       data linear in t, smooth in space, gap h^2, h = 1/4, 1/8, k = 1, 2.
%
% Pieces: set the environment variable HDG_P2_PARTS to a subset of
% "abcde" (default all), HDG_P2_DEGREES to the degrees of (b) (e.g. "1 2";
% default 1:3) and HDG_P2_VARIANTS to the variants of (c) (e.g. "euler";
% default "euler bdf2 bdf2_exact_start"), e.g.
%   setenv("HDG_P2_PARTS", "c"); setenv("HDG_P2_VARIANTS", "euler");
%   test_stokes_biot_transient
clc; clear;
parts = string(getenv("HDG_P2_PARTS"));
if parts == "", parts = "abcde"; end
degrees_b = str2num(getenv("HDG_P2_DEGREES")); %#ok<ST2NM>
if isempty(degrees_b), degrees_b = 1:3; end
variants_c = split(strtrim(string(getenv("HDG_P2_VARIANTS")))).';
if variants_c == "", variants_c = ["euler", "bdf2", "bdf2_exact_start"]; end
names = ["u_F", "sigma_F", "eta", "v", "sigma_P", "q", "p"];
base = {"storage", 1.0, "friction_bjs", 1.0};
n_fail = 0;

%% (a) bookkeeping
if contains(parts, "a")
    F = struct("velocity", 3*ones(2), "velocity_previous", ones(2));
    [cdt, h] = fluid_step_history("bdf2", 0.3, F);
    ok = abs(cdt - 0.2) < 1e-15 && max(abs(h.velocity - (11/3)*ones(2)), [], "all") < 1e-14 ...
        && isequal(h.convective_velocity, {});
    n_fail = n_fail + check(ok, "fluid_step_history bdf2: cdt = 2 dt/3, u* = (4 u^n - u^{n-1})/3");
    [cdt, h] = fluid_step_history("euler", 0.3, F);
    n_fail = n_fail + check(cdt == 0.3 && isequal(h.velocity, 3*ones(2)), ...
        "fluid_step_history euler: cdt = dt, u* = u^n");
    [cdt, h] = fluid_step_history("bdf2", 0.3, struct("velocity", 3*ones(2), "velocity_previous", []));
    n_fail = n_fail + check(cdt == 0.3 && isequal(h.velocity, 3*ones(2)), ...
        "fluid_step_history bdf2 without level n-1 ([]): the Euler step");
    [cdt, h] = fluid_step_history("bdf2", 0.3, struct("velocity", 3*ones(2)));
    n_fail = n_fail + check(cdt == 0.3 && isequal(h.velocity, 3*ones(2)), ...
        "fluid_step_history bdf2 without the field velocity_previous: the Euler step");
    [cdt, h] = fluid_step_history("bdf2", 0.3, struct("velocity", [], "velocity_previous", 0));
    n_fail = n_fail + check(abs(cdt - 0.2) < 1e-15 && isequal(h.velocity, 0), ...
        "fluid_step_history bdf2, zero levels: u* = the scalar 0");
    [cdt, h] = fluid_step_history("bdf2", 0.3);
    n_fail = n_fail + check(cdt == 0.3 && isequal(h.velocity, 0) && isequal(h.convective_velocity, {}), ...
        "fluid_step_history, empty state: Euler, zero history, no convective velocity");
    U = struct("velocity", ones(2), "velocity_previous", [], ...
        "u_star", {{2*ones(3); 4*ones(3)}}, "u_star_previous", {{ones(3); 3*ones(3)}});
    [~, h] = fluid_step_history("bdf2", 0.3, U);
    ok = isequal(h.convective_velocity, {3*ones(3); 5*ones(3)});
    [~, h] = fluid_step_history("euler", 0.3, U);
    ok = ok && isequal(h.convective_velocity, {2*ones(3); 4*ones(3)});
    U.u_star_previous = [];
    [~, h] = fluid_step_history("bdf2", 0.3, U);
    ok = ok && isequal(h.convective_velocity, {2*ones(3); 4*ones(3)});
    n_fail = n_fail + check(ok, "fluid_step_history convective_velocity: 2u*^n - u*^{n-1} (bdf2), u*^n (euler, no level n-1)");

    % initial_state and advance_state on one coupled step, k = 1, h = 1/4
    k = 1; h = 1/4; dt = 0.1; t0 = 0.5;
    f = patch_fields(k, "linear");
    mms = StokesBiotCoupling.get_manufactured_solution(["0"; "1"], f{:}, base{:}, ...
        density_fluid = 1, density_biot = 1);
    [fm, bm, curve] = flat_meshes(h, 0);
    for scheme = ["euler", "bdf2"]
        [fp, bp] = make_problems(k, mms, dt, scheme, t0);
        s0 = StokesBiotCoupling.initial_state(fp, bp, fm, bm, Order = k);
        d0 = mms.at(t0);
        biot_only = initial_state(bp, bm, k);
        ok = s0.time == t0 && s0.step == 0 && isempty(s0.fluid.velocity_previous) ...
            && max(abs(s0.fluid.velocity - L2proj2D(fm, d0.u_F, k)), [], "all") < 1e-14 ...
            && isequal(s0.biot, biot_only.biot);
        n_fail = n_fail + check(ok, scheme + ": StokesBiotCoupling.initial_state (Biot state + projected u^0)");
        s0.fluid.u_star = {ones(3, 2); ones(3, 2)};
        [fp, bp] = set_data(fp, bp, mms.at(t0 + dt));
        [fs, bs, s1] = StokesBiotCoupling.solve(fp, bp, fm, bm, curve, mms.at(t0 + dt), ...
            Order = k, friction_bjs = 1, state = s0);
        ok = abs(s1.time - (t0 + dt)) < 1e-15 && s1.step == 1 ...
            && isequal(s1.fluid.velocity, fs.Unknowns.field(2).dof) ...
            && isequal(s1.fluid.velocity_previous, s0.fluid.velocity) ...
            && isequal(s1.fluid.velocity_skeleton, fs.Unknowns.hybrid(1).dof) ...
            && isequal(s1.fluid.stress, fs.Unknowns.field(1).dof) ...
            && isequal(size(s1.fluid.pressure), [3, size(fm.elements, 2)]) ...
            && isequal(s1.fluid.u_star, s0.fluid.u_star) ...
            && isequal(s1.biot.displacement, bs.Unknowns.field(2).dof);
        n_fail = n_fail + check(ok, scheme + ": advance_state (time, step once; fluid fields; u_star kept)");
        % the fluid pressure: -tr(sigma)/2 of a degree-k patch is exact
        d1 = mms.at(t0 + dt);
        e = compute_l2_error(scalar_field(fs.Unknowns.field(3), s1.fluid.pressure), fm, ...
            Datum(d1.p_F));
        n_fail = n_fail + check(e < 1e-10, sprintf("%s: state.fluid.pressure = p_F (error %.1e)", scheme, e));
    end
    % one step for both media
    [fp, bp] = make_problems(k, mms, dt, "bdf2", t0);
    fp.Data.time_step = 2*dt;
    n_fail = n_fail + check(throws(@() StokesBiotCoupling.interface_problem(fp, bp, fm, bm, ...
        curve, mms.at(t0 + dt), Order = k), "StokesBiotCoupling:timeStep"), ...
        "rho_F > 0: a fluid time_step other than the Biot one is an error");
    fp.Data.time_step = dt;
    s = s1;                                % has both levels n - 1
    s.fluid.velocity_previous = [];
    n_fail = n_fail + check(throws(@() StokesBiotCoupling.interface_problem(fp, bp, fm, bm, ...
        curve, mms.at(t0 + 2*dt), Order = k, state = s), "StokesBiotCoupling:historyLevels"), ...
        "bdf2: only the Biot state has the level n-1: error historyLevels");
    s = s1;
    s.biot.displacement_skeleton_previous = [];
    n_fail = n_fail + check(throws(@() StokesBiotCoupling.interface_problem(fp, bp, fm, bm, ...
        curve, mms.at(t0 + 2*dt), Order = k, state = s), "StokesBiotCoupling:historyLevels"), ...
        "bdf2: only the fluid state has the level n-1: error historyLevels");

    % V0: rho_F = 0 does not read state.fluid ...
    mms0 = StokesBiotCoupling.get_manufactured_solution(["0"; "1"], f{:}, base{:});
    [fp, bp] = make_problems(k, mms0, dt, "bdf2", t0);
    fp.Data.density = 0;
    s = StokesBiotCoupling.initial_state(fp, bp, fm, bm, Order = k);
    [fp, bp] = set_data(fp, bp, mms0.at(t0 + dt));
    [fs_a, bs_a] = StokesBiotCoupling.solve(fp, bp, fm, bm, curve, mms0.at(t0 + dt), ...
        Order = k, friction_bjs = 1, state = s);
    s.fluid.velocity = rand(size(s.fluid.velocity));
    s.fluid.velocity_previous = rand(size(s.fluid.velocity));
    [fs_b, bs_b] = StokesBiotCoupling.solve(fp, bp, fm, bm, curve, mms0.at(t0 + dt), ...
        Order = k, friction_bjs = 1, state = s);
    ok = isequal(fs_a.Unknowns.field(2).dof, fs_b.Unknowns.field(2).dof) ...
        && isequal(bs_a.Unknowns.field(2).dof, bs_b.Unknowns.field(2).dof);
    n_fail = n_fail + check(ok, "rho_F = 0: the step does not read state.fluid");
    % ... and reproduces the code before the transient branch (errors of
    % stokes_biot_case recorded there: flat, k = 2, h = 1/4, gap 0, one
    % step from rest; then bdf2, exact rate, 3 steps, dt = 0.1, gap h^2)
    recorded = {[0.0127753488298371 0.00454697023976069 0.00221172423538017 ...
        0.000562040837357954 0.00689622095438911 0.00311170383969279], [0.0433595658531354 0.00486940814903846 0.0023724966233459 ...
        0.000535921547678854 0.0062222425088327 0.00213619170990052]};
    runs = {stokes_biot_case("flat", 2, 1/4, gap = 0), ...
        stokes_biot_case("flat", 2, 1/4, time_scheme = "bdf2", rate = "exact", ...
            n_steps = 3, time_step = 0.1, storage = 1)};
    for i = 1:2
        e = cell2mat(struct2cell(runs{i}.errors)).';
        dev = max(abs(e - recorded{i})./recorded{i});
        n_fail = n_fail + check(dev <= 1e-10, sprintf( ...
            "V0 case %d: rho_F = 0 errors as before the transient branch (max rel. deviation %.1e)", i, dev));
    end
end

%% (b) patch test
if contains(parts, "b")
    tol = 1e-10;
    cases = struct("label", {"euler, linear in t", "bdf2 (Euler start), linear in t", ...
        "bdf2 (exact start), quadratic in t"}, "scheme", {"euler", "bdf2", "bdf2"}, ...
        "time", {"linear", "linear", "quadratic"}, "start", {"initial", "initial", "exact"}, ...
        "n_steps", {3, 3, 2});
    fprintf("\n(b) patch test, flat Sigma, no gap, rho_F = rho_P = 1, dt = 0.1, from t = 1/2\n");
    for k = degrees_b
        for c = cases
            f = patch_fields(k, c.time);
            mms = StokesBiotCoupling.get_manufactured_solution(["0"; "1"], ...
                f{:}, base{:}, density_fluid = 1, density_biot = 1);
            for h = [1/4, 1/6]
                out = march(k, h, 0, mms, c.scheme, 0.1, 0.5, c.n_steps, c.start);
                [worst, at] = max(out.relative(end, :));
                fprintf("  k %d  h 1/%d  %-34s max rel. error %.2e (%s)\n", k, round(1/h), ...
                    c.label, worst, names(at));
                n_fail = n_fail + check(worst <= tol, sprintf("(b) k = %d, h = 1/%d, %s: %s %.2e", ...
                    k, round(1/h), c.label, names(at), worst));
            end
        end
    end
end

%% (c) time order
if contains(parts, "c")
    k = 2; h = 1/4; T = 1;
    dts = 1./[4, 8, 16, 32, 64, 128];
    thresholds = struct("euler", 0.9, "bdf2", -Inf, "bdf2_exact_start", 1.8);
    mms = StokesBiotCoupling.get_manufactured_solution(["0"; "1"], base{:}, ...
        stream_function = "(2 + cos(2*t))*(x^3 - 2*x*y^2 + x^2*y + y^3/3 + x*y)", ...
        pressure_fluid = "(1 + sin(t))*(x^2 - x*y + 2*y^2 + x)", ...
        displacement = ["exp(t)*(x^2 + x*y + y)"; "(2 - cos(3*t))*(x*y - y^2 + x)"], ...
        pressure_biot = "cos(2*t)*(x^2 - x*y + y^2 + 1)", ...
        density_fluid = 1, density_biot = 1);
    for variant = variants_c
        scheme = extractBefore(variant + "_", "_");
        start = "initial";
        if variant == "bdf2_exact_start", start = "exact"; end
        errors = zeros(numel(dts), numel(names));
        for i = 1:numel(dts)
            n_steps = round(T/dts(i));
            if start == "exact", n_steps = n_steps - 1; end
            out = march(k, h, h^2, mms, scheme, dts(i), 0, n_steps, start);
            errors(i, :) = out.absolute(end, :);
        end
        eoc = log2(errors(1:end-1, :)./errors(2:end, :));
        fprintf("\n(c) %s: errors at T = %g, k = %d, h = 1/%d, gap h^2\n", ...
            strrep(variant, "_", " "), T, k, round(1/h));
        print_table("dt", dts, errors, eoc, names);
        threshold = thresholds.(variant);
        if isfinite(threshold)
            n_fail = n_fail + check(all(eoc(end, :) >= threshold), ...
                sprintf("(c) %s: EOC in time of all fields (two finest steps) >= %.1f", ...
                strrep(variant, "_", " "), threshold));
        else
            fprintf("  (reported) %s: EOC of the two finest steps %s\n", strrep(variant, "_", " "), ...
                strjoin(names + " " + compose("%.2f", eoc(end, :)), ", "));
        end
    end
end

%% (d) V5: three BDF2 steps
if contains(parts, "d")
    k = 2; h = 1/8; dt = 1/20; t0 = 0.5;
    mms = StokesBiotCoupling.get_manufactured_solution(["0"; "1"], base{:}, ...
        stream_function = "(2 + t - t^2)*sin(pi*x)*sin(pi*y)", ...
        pressure_fluid = "(1 + t^2)*sin(pi*x)*exp(y - 0.5)", ...
        displacement = ["(1 + t + t^2/2)*x*sin(y)"; "(2 - t + 2*t^2)*y*cos(x)"], ...
        pressure_biot = "(1 - 2*t^2)*cos(pi*x)*sin(pi*y)", ...
        density_fluid = 1, density_biot = 1);
    for start = ["exact", "initial"]
        out = march(k, h, h^2, mms, "bdf2", dt, t0, 3, start);
        ratio = out.absolute(end, :)./out.absolute(1, :);
        fprintf("\n(d) bdf2, %s start, k = %d, h = 1/%d, gap h^2, dt = %g, from t = %g\n", ...
            start, k, round(1/h), dt, t0);
        print_table("step", 1:3, out.absolute, [], names);
        fprintf("  error(step 3)/error(step 1): %s\n", strjoin(names + " " + compose("%.2f", ratio), ", "));
        if start == "exact"
            n_fail = n_fail + check(all(ratio <= 2), ...
                "(d) V5: every field's error after 3 bdf2 steps within 2x the first step's");
        end
    end
end

%% (e) sanity space rate
if contains(parts, "e")
    t0 = 0.5; dt = 0.25;
    mms = StokesBiotCoupling.get_manufactured_solution(["0"; "1"], base{:}, ...
        stream_function = "(1 + t)*sin(pi*x)*sin(pi*y)", ...
        pressure_fluid = "(2 - t)*sin(pi*x)*exp(y - 0.5)", ...
        displacement = ["(1 + t)*x*sin(y)"; "(2 - t)*y*cos(x)"], ...
        pressure_biot = "(1 + 2*t)*cos(pi*x)*sin(pi*y)", ...
        density_fluid = 1, density_biot = 1);
    hs = [1/4, 1/8];
    for k = 1:2
        errors = zeros(numel(hs), numel(names));
        for i = 1:numel(hs)
            out = march(k, hs(i), hs(i)^2, mms, "euler", dt, t0, 1, "initial");
            errors(i, :) = out.absolute(end, :);
        end
        fprintf("\n(e) k = %d, one Euler step from the exact state (linear in t), gap h^2 (reported)\n", k);
        print_table("h", hs, errors, log2(errors(1, :)./errors(2, :)), names);
    end
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

function f = patch_fields(k, time)
% polynomial fields of degree <= k in space (stream function k + 1), with
% a time factor linear or quadratic in t, as name-value pairs of
% StokesBiotCoupling.get_manufactured_solution
switch time
    case "linear",    g = ["(1 + 2*t)", "(2 - t)", "(1 + t)", "(3 - 2*t)", "(1 + 3*t)"];
    case "quadratic", g = ["(1 + 2*t + t^2)", "(2 - t + 2*t^2)", "(1 + t - t^2)", ...
                           "(3 - 2*t + t^2)", "(1 + 3*t + 2*t^2)"];
end
f = {"stream_function", g(1) + "*(" + poly_string(k + 1, 1) + ")", ...
     "pressure_fluid",  g(2) + "*(" + poly_string(k, 2) + ")", ...
     "displacement",   [g(3) + "*(" + poly_string(k, 3) + ")"; g(4) + "*(" + poly_string(k, 4) + ")"], ...
     "pressure_biot",   g(5) + "*(" + poly_string(k, 5) + ")"};
end

function s = poly_string(degree, seed)
% sum of c_ij x^i y^j over i + j <= degree, fixed non-trivial coefficients
terms = strings(0);
for i = 0:degree
    for j = 0:degree - i
        c = (-1)^(i + seed) * (1 + mod(3*i + 5*j + seed, 7)) / 4;
        terms(end + 1) = sprintf("(%g)*x^%d*y^%d", c, i, j); %#ok<AGROW>
    end
end
s = join(terms, " + ");
end

function [fluid_mesh, biot_mesh, curve] = flat_meshes(h, gap)
curve = FlatCurveHorizontal(0.5);
[fluid_mesh, biot_mesh] = fitted_split(curve, [0.8*h, h], gap);
end

function [fp, bp] = make_problems(k, mms, dt, scheme, t)
% the fluid and Biot problems with the initial data at time t
P = mms.parameters;
d = mms.at(t);
fp = IncompressibleElasticityProblem(k, viscosity = P.viscosity, ...
    density = P.density_fluid, time_step = dt, source_term = d.f_F, ...
    datum_dirichlet = d.u_F, pressure_mean = "pinned", initial_velocity = d.u_F);
bp = BiotProblem(k, young_modulus = P.young_modulus, poisson_ratio = P.poisson_ratio, ...
    biot_alpha = P.biot_alpha, storage = P.storage, diffusion = P.diffusion, ...
    hooke_modulus = P.hooke_modulus, density = P.density_biot, ...
    time = t, time_step = dt, time_scheme = scheme, ...
    source_term = d.f_P, source_mass = d.g_P, initial_displacement = d.u_P, ...
    initial_velocity = d.u_t_P, initial_pressure = d.p_P, ...
    datum_dirichlet = {d.u_P{1}; d.u_P{2}; d.p_P});
end

function [fp, bp] = set_data(fp, bp, d)
% sources and Dirichlet data at t_{n+1}
fp.Data.source_term = d.f_F;
fp.Data.datum_dirichlet = d.u_F;
bp.Data.source_term = d.f_P;
bp.Data.source_mass = d.g_P;
bp.Data.datum_dirichlet = {d.u_P{1}; d.u_P{2}; d.p_P};
end

function out = march(k, h, gap, mms, scheme, dt, t0, n_steps, start)
% n_steps coupled steps from t0 (start = "initial": the projected state
% at t0, a bdf2 march starting with its Euler step; "exact": the
% projected states at t0 and t0 + dt as levels n-1 and n, the march from
% t0 + dt). out.absolute / out.relative: errors of the fields after each
% step (rows), relative to the L2 norm of the exact field
[fm, bm, curve] = flat_meshes(h, gap);
[fp, bp] = make_problems(k, mms, dt, scheme, t0);
state = StokesBiotCoupling.initial_state(fp, bp, fm, bm, Order = k);
if start == "exact"
    previous = state;
    [fp, bp] = make_problems(k, mms, dt, scheme, t0 + dt);
    state = StokesBiotCoupling.initial_state(fp, bp, fm, bm, Order = k);
    for name = ["displacement_skeleton", "stress", "rotation", "pressure", ...
            "displacement", "velocity"]
        state.biot.(name + "_previous") = previous.biot.(name);
    end
    state.fluid.velocity_previous = previous.fluid.velocity;
    state.step = 1;
end
[out.absolute, out.relative] = deal(zeros(n_steps, 7));
for n = 1:n_steps
    d = mms.at(state.time + dt);
    [fp, bp] = set_data(fp, bp, d);
    [fs, bs, state] = StokesBiotCoupling.solve(fp, bp, fm, bm, curve, d, ...
        Order = k, friction_bjs = mms.parameters.friction_bjs, state = state);
    [out.absolute(n, :), out.relative(n, :)] = field_errors(fs, bs, state, fm, bm, mms.at(state.time));
end
end

function [e, rel] = field_errors(fs, bs, state, fm, bm, d)
% L2 errors of u_F, sigma_F, eta, v, sigma_P, q, p and the same relative
% to the L2 norm of the exact field (absolute where it is below 1e-12)
ff = fs.Unknowns.field;
bf = bs.Unknowns.field;
v = bf(2);
v.dof = state.biot.velocity;
list = {ff(2), fm, d.u_F; ff(1), fm, d.sigma_F; bf(2), bm, d.u_P; v, bm, d.u_t_P; ...
        bf(1), bm, d.sigma_P; bf(4), bm, d.q_P; bf(5), bm, d.p_P};
[e, rel] = deal(zeros(1, size(list, 1)));
for i = 1:size(list, 1)
    [u, mesh, exact] = list{i, :};
    e(i) = compute_l2_error(u, mesh, Datum(exact));
    z = u;
    z.dof = 0*z.dof;
    rel(i) = e(i)/max(compute_l2_error(z, mesh, Datum(exact)), 1e-12);
end
end

function p = scalar_field(template, dof)
% the scalar Dubiner field of the dofs dof, on the mesh and space of
% template (the fluid rotation, of the same degree)
p = template;
p.dof = dof;
end

function print_table(label, steps, errors, eoc, names)
fprintf("%8s", label);
fprintf(" %10s %5s", [names; repmat("EOC", 1, numel(names))]);
fprintf("\n");
for i = 1:numel(steps)
    fprintf("%8.4g", steps(i));
    for j = 1:numel(names)
        if i == 1 || isempty(eoc)
            fprintf(" %10.3e %5s", errors(i, j), "-");
        else
            fprintf(" %10.3e %5.2f", errors(i, j), eoc(i-1, j));
        end
    end
    fprintf("\n");
end
end
