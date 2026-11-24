% test_stokes_biot_transient_patch.m
% Transient Stokes-Biot coupling (StokesBiotCoupling.solve, rho_F = rho_P
% = 1, no convection) beyond the flat, gap-free patch test of
% test_stokes_biot_transient: gaps, the circle, long marches, and the time
% order on the circle. Unit square, all-Dirichlet outer boundaries,
% pinned fluid pressure, s0 = 1, friction_bjs = 1, alpha = 1, kappa = I.
%   flat:   Sigma y = 1/2, fluid above, fitted_split (0.8 h, h), gap delta
%   circle: centre (1/2, 1/2), R = 1/4, Biot inside, fluid outside,
%           interface_meshes(curve, "fitted", 0.8 h, gap = delta,
%           h_solid = h) (both meshes are polygons; neither reaches Sigma)
%
%   (a) patch test, k = 1, 2, 3: every field a polynomial of degree <= k
%       in space (stream function k + 1), times a factor in t:
%         euler, linear in t, 3 steps from the initial state at t = 1/2;
%         bdf2 from an exact two-level start, quadratic in t, 2 steps;
%       cases flat gap h^2 and gap h (h = 1/4, 1/8), circle gap 0 and
%       h^2 (h = 1/8, 1/12). Every field after every step is checked to
%       <= 1e-10 relative (to the L2 norm of the exact field), except gap
%       h at k = 3, which is reported: there the error is set at the
%       first step by the conditioning of the Taylor extrapolation over a
%       path of length delta = h (1.7e-10 for q at h = 1/4, delta = 1/4;
%       the quasi-static test_stokes_biot_patch_gap gives 9.4e-11 there)
%       and stays at that level over the steps (see (b) for long
%       marches). The error after every step is printed for gap h.
%   (b) long marches, reported: flat gap h, patch fields, k = 2, 3, h =
%       1/8, 1/16: 20 euler steps (linear in t) and 20 bdf2 steps from
%       the exact start (quadratic in t), dt = 0.1. The max relative
%       error after every step and its growth factor per step.
%   (c) time order on the circle, gap h^2, k = 2, h = 1/8: stream
%       function of degree 3 and p_F, eta, p_P of degree 2 in space
%       (in the discrete space: the error is the time error up to the
%       patch-test level), smooth in t, T = 1, dt = 1/8 ... 1/128,
%       variants euler, bdf2 (the toolbox's Euler start) and bdf2 from an
%       exact two-level start. EOC of all fields between the two finest
%       steps: >= 0.9 (euler), >= 1.8 (bdf2 exact start), reported (bdf2
%       Euler start).
%
% Pieces: environment variables HDG_TP_PARTS (subset of "abc", default
% all), HDG_TP_DEGREES (degrees of (a) and (b), e.g. "1 2"; default 1:3
% for (a), 2:3 for (b)), HDG_TP_CASES (cases of (a), default
% "flat_h2 flat_h circle_zero circle_h2") and HDG_TP_VARIANTS (variants of
% (c), default "euler bdf2 bdf2_exact_start"), e.g.
%   setenv("HDG_TP_PARTS", "a"); setenv("HDG_TP_DEGREES", "3");
%   test_stokes_biot_transient_patch
clc; clear;
parts = string(getenv("HDG_TP_PARTS"));
if parts == "", parts = "abc"; end
degrees = str2num(getenv("HDG_TP_DEGREES")); %#ok<ST2NM>
cases_a = split(strtrim(string(getenv("HDG_TP_CASES")))).';
if cases_a == "", cases_a = ["flat_h2", "flat_h", "circle_zero", "circle_h2"]; end
variants_c = split(strtrim(string(getenv("HDG_TP_VARIANTS")))).';
if variants_c == "", variants_c = ["euler", "bdf2", "bdf2_exact_start"]; end
names = ["u_F", "sigma_F", "eta", "v", "sigma_P", "q", "p"];
asserted_cases = ["flat_h2", "flat_h", "circle_zero", "circle_h2"];
tol = 1e-10;
n_fail = 0;

%% (a) patch test
if contains(parts, "a")
    degrees_a = degrees;
    if isempty(degrees_a), degrees_a = 1:3; end
    schemes = struct("label", {"euler, linear in t", "bdf2 (exact start), quadratic in t"}, ...
        "scheme", {"euler", "bdf2"}, "time", {"linear", "quadratic"}, ...
        "start", {"initial", "exact"}, "n_steps", {3, 2});
    fprintf("\n(a) patch test, rho_F = rho_P = 1, dt = 0.1, from t = 1/2\n");
    for k = degrees_a
        for c = cases_a
            [problem, gap_rule] = split_case(c);
            hs = [1/4, 1/8];
            if problem == "circle", hs = [1/8, 1/12]; end
            for s = schemes
                f = patch_fields(k, s.time);
                mms = StokesBiotCoupling.get_manufactured_solution(normal_of(problem), ...
                    f{:}, storage = 1, friction_bjs = 1, ...
                    density_fluid = 1, density_biot = 1);
                for h = hs
                    out = march(k, problem, h, gap_of(gap_rule, h), mms, s.scheme, 0.1, 0.5, ...
                        s.n_steps, s.start);
                    [worst, at] = max(out.relative, [], 2);
                    fprintf("  k %d  %-11s h 1/%-2d %-36s max rel. error %.2e (%s)\n", k, c, ...
                        round(1/h), s.label, worst(end), names(at(end)));
                    if gap_rule == "h"
                        for n = 1:s.n_steps
                            fprintf("      step %d: max rel. %.2e (%s)%s\n", n, worst(n), ...
                                names(at(n)), sprintf(" %s %.1e", [names; out.relative(n, :)]));
                        end
                    end
                    label = sprintf("(a) k = %d, %s, h = 1/%d, %s: %s %.2e", k, c, ...
                        round(1/h), s.label, names(at(end)), max(worst));
                    if ismember(c, asserted_cases) && ~(c == "flat_h" && k == 3)
                        n_fail = n_fail + check(max(worst) <= tol, label);
                    else
                        fprintf("  (reported) %s\n", label);
                    end
                end
            end
        end
    end
end

%% (b) long marches at gap h
if contains(parts, "b")
    degrees_b = degrees;
    if isempty(degrees_b), degrees_b = 2:3; end
    n_steps = 20;
    fprintf("\n(b) long marches, flat gap h, patch fields, dt = 0.1, %d steps (reported)\n", n_steps);
    for k = degrees_b
        for s = struct("label", {"euler, linear in t", "bdf2 (exact start), quadratic in t"}, ...
                "scheme", {"euler", "bdf2"}, "time", {"linear", "quadratic"}, ...
                "start", {"initial", "exact"})
            f = patch_fields(k, s.time);
            mms = StokesBiotCoupling.get_manufactured_solution(["0"; "1"], ...
                f{:}, storage = 1, friction_bjs = 1, ...
                density_fluid = 1, density_biot = 1);
            for h = [1/8, 1/16]
                out = march(k, "flat", h, h, mms, s.scheme, 0.1, 0.5, n_steps, s.start);
                [worst, at] = max(out.relative, [], 2);
                growth = worst(2:end)./worst(1:end-1);
                fprintf("  k %d  h 1/%-2d %s: max rel. error per step (field), growth factor\n", ...
                    k, round(1/h), s.label);
                for n = 1:n_steps
                    g = NaN;
                    if n > 1, g = growth(n - 1); end
                    fprintf("    %2d  %.3e (%-7s) %6.2f\n", n, worst(n), names(at(n)), g);
                end
                fprintf("  (reported) k %d h 1/%d %s: step 1 %.2e, step %d %.2e, geometric mean growth %.2f\n", ...
                    k, round(1/h), s.label, worst(1), n_steps, worst(end), ...
                    (worst(end)/worst(1))^(1/(n_steps - 1)));
            end
        end
    end
end

%% (c) time order on the circle
if contains(parts, "c")
    k = 2; h = 1/8; T = 1;
    dts = 1./[8, 16, 32, 64, 128];
    thresholds = struct("euler", 0.9, "bdf2", -Inf, "bdf2_exact_start", 1.8);
    mms = StokesBiotCoupling.get_manufactured_solution(normal_of("circle"), ...
        storage = 1, friction_bjs = 1, ...
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
            out = march(k, "circle", h, h^2, mms, scheme, dts(i), 0, n_steps, start);
            errors(i, :) = out.absolute(end, :);
        end
        eoc = log2(errors(1:end-1, :)./errors(2:end, :));
        fprintf("\n(c) circle, %s: errors at T = %g, k = %d, h = 1/%d, gap h^2\n", ...
            strrep(variant, "_", " "), T, k, round(1/h));
        print_table("dt", dts, errors, eoc, names);
        threshold = thresholds.(variant);
        if isfinite(threshold)
            n_fail = n_fail + check(all(eoc(end, :) >= threshold), ...
                sprintf("(c) circle %s: EOC in time of all fields (two finest steps) >= %.1f", ...
                strrep(variant, "_", " "), threshold));
        else
            fprintf("  (reported) circle %s: EOC of the two finest steps %s\n", ...
                strrep(variant, "_", " "), strjoin(names + " " + compose("%.2f", eoc(end, :)), ", "));
        end
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

function [problem, gap_rule] = split_case(c)
parts = split(c, "_");
problem = parts(1);
gap_rule = parts(2);
end

function delta = gap_of(rule, h)
switch rule
    case "zero", delta = 0;
    case "h2",   delta = h^2;
    case "h",    delta = h;
end
end

function n = normal_of(problem)
% nu_p, outward from the Biot domain
if problem == "flat"
    n = ["0"; "1"];
else
    r = "sqrt((x - 0.5)^2 + (y - 0.5)^2)";
    n = ["(x - 0.5)/" + r; "(y - 0.5)/" + r];
end
end

function f = patch_fields(k, time)
% polynomial fields of degree <= k in space (stream function k + 1), with
% a time factor linear or quadratic in t (as test_stokes_biot_transient)
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

function [fluid_mesh, biot_mesh, curve] = meshes(problem, h, gap)
if problem == "flat"
    curve = FlatCurveHorizontal(0.5);
    [fluid_mesh, biot_mesh] = fitted_split(curve, [0.8*h, h], gap);
else
    curve = CircularCurve([0.5; 0.5], 0.25);
    [fluid_mesh, biot_mesh] = interface_meshes(curve, "fitted", 0.8*h, ...
        box = [0 1 0 1], fluid_side = "outer", gap = gap, h_solid = h);
end
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

function out = march(k, problem, h, gap, mms, scheme, dt, t0, n_steps, start)
% n_steps coupled steps from t0 (start = "initial": the projected state
% at t0, a bdf2 march starting with its Euler step; "exact": the
% projected states at t0 and t0 + dt as levels n-1 and n, the march from
% t0 + dt). out.absolute / out.relative: errors of the fields after each
% step (rows), relative to the L2 norm of the exact field
[fm, bm, curve] = meshes(problem, h, gap);
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
    fp.Data.source_term = d.f_F;
    fp.Data.datum_dirichlet = d.u_F;
    bp.Data.source_term = d.f_P;
    bp.Data.source_mass = d.g_P;
    bp.Data.datum_dirichlet = {d.u_P{1}; d.u_P{2}; d.p_P};
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
