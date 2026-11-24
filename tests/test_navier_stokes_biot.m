% test_navier_stokes_biot.m
% TPM-coupled Navier-Stokes-Biot step (StokesBiotCoupling.solve_navier_stokes):
% fluid and Biot inertia (rho_F = rho_P = 1), the fluid's convective term
% rho_F (beta . grad) u with beta = u* (divergence_free_velocity) of a
% fluid solution, modes "iterated" (Picard) and "semi-implicit" (one solve
% with the predictor u*^n, or 2 u*^n - u*^{n-1} under bdf2). Manufactured
% data of StokesBiotCoupling.get_manufactured_solution(..., density_fluid
% = 1, density_biot = 1, convection = true); s0 = 1, friction_bjs = 1,
% Dirichlet data on the whole outer boundary, pinned fluid pressure, mu =
% 1. The fluid mesh.tau is set per face to 1 + rho B_f / mu with B_f the
% max of |u_F| over the face (Gauss points) and over the time window, so
% S - rho beta . n/2 >= mu + rho (B_f - beta . n/2) > 0 for beta = u* ~ u_F
% (info.s_margin is checked > 0 in every run).
% Flat Sigma: y = 1/2, fitted_split (0.8 h, h) with the gap; circle: the
% circle of radius 1/4 about (1/2, 1/2), interface_meshes "fitted" (0.8 h,
% h_solid = h) with gap h^2.
%
%   (a) bookkeeping, k = 2, flat h = 1/4: semi-implicit with a state that
%       has u_star = the dof-identical manual interface_problem + solve
%       with beta = that u_star (euler: u*^n; bdf2 with both levels: 2u*^n
%       - u*^{n-1}); the returned state (time, step once; u_star,
%       u_star_previous; fluid fields); the first step without u_star (a
%       Stokes pre-solve, u_star_previous = []); the Cache (iterate and
%       step) gives identical results; the error
%       StokesBiotCoupling:stabilization of interface_problem; s_margin,
%       max_bn, max_bn_interface against an independent evaluation with
%       the exact (polynomial) u_F at the same points; the warning
%       StokesBiotCoupling:stabilization with a too small tau.
%   (b) patch test, Euler, linear in t, 3 steps from the initial state at
%       t = 1/2, dt = 0.1, iterated to 1e-12: flat gap 0 and gap h^2 (h =
%       1/4, 1/6), k = 1, 2, 3, every field <= 1e-10 relative; circle gap
%       h^2 (h = 1/8, 1/12) reported. For k = 3 the stream function has
%       degree 3 (u_F of degree 2), the other fields degree 3: the
%       convective source (u . grad) u of degree 2 deg(u) - 1 times a
%       test function of degree k must stay within the element source
%       rule (exact to degree 2k + 1 = 7), so deg(u) <= 2 at k = 3.
%   (c) V5: 3 bdf2 steps from an exact two-level start (projections of
%       the exact fields at t0, t0 + dt; u_star, u_star_previous = their
%       P_{k+1} projections), data smooth in space and quadratic in t, k =
%       2, flat h = 1/8, gap h^2, dt = 1/20, iterated to 1e-10: every
%       field's error after 3 steps within 2x the error after step 1. The
%       Picard-converged bdf2 step is exact in time for quadratic-in-t
%       data: the exact solution satisfies the discrete-in-time equation
%       (bdf2 differentiates quadratics exactly, and the convective term
%       and source are both taken at t_{n+1}), so only the space error
%       remains.
%   (d) V6 time order, k = 2, flat h = 1/4, gap h^2, smooth in t, T = 1,
%       dt = 1/4 ... 1/64: EOC of the two finest steps, all fields: euler
%       >= 0.9 (iterated and semi-implicit), bdf2 from an exact start >=
%       1.8 (iterated and semi-implicit), bdf2 from the toolbox's Euler
%       start reported (both modes).
%   (e) Picard vs semi-implicit: from the runs of (d), the L2 difference
%       of the two modes' fields at T per dt and its rate; Picard
%       iterations per step (reported).
%   (f) V7 space: 2 bdf2 steps from an exact start, data quadratic in t
%       (no time error, see (c)), dt = 0.2, iterated to 1e-10: flat gap
%       h^2 (h = 1/4, 1/8, 1/16) and circle gap h^2 (h = 1/10, 1/20,
%       1/40); rates coarsest-finest >= k + 1 - 0.3 (u_F, eta, p), k + 1 -
%       0.5 (sigma_F, v, sigma_P, q) flat, - 0.4 / - 0.6 circle. The step
%       is not small on purpose: the start state is the projection of the
%       exact fields, not a discrete solution, and their O(h^(k+1))
%       difference enters v (and eta through it) amplified by 1/cdt (the
%       start transient of test_stokes_biot_transient (d)). With dt =
%       0.05 the flat ladder gives v 1.26 (k = 1), eta 2.56 / v 2.37 (k =
%       2), eta 3.51 (k = 3), the same with and without convection; with
%       dt = 0.2 the rates are optimal.
%
% Pieces: environment variables NSB_PARTS (subset of "abcdef", default
% all), NSB_DEGREES (degrees of (b) and (f), e.g. "1,2"; default 1:3 for
% (b), 1:2 for (f)), NSB_VARIANTS (variants of (d), default "euler bdf2
% bdf2_exact_start"), NSB_MODES (modes of (d), default "iterated
% semi-implicit"), NSB_LADDERS (ladders of (f), default "flat circle"),
% e.g.
%   setenv("NSB_PARTS", "f"); setenv("NSB_DEGREES", "3"); test_navier_stokes_biot
% (e) needs (d) in the same run. k = 0 is skipped (no curl bubble in the
% stress elements).
clc; clear;
addpath(fullfile(fileparts(mfilename("fullpath")), "helpers"));
parts = string(getenv("NSB_PARTS"));
if parts == "", parts = "abcdef"; end
degrees = str2double(split(string(getenv("NSB_DEGREES")), ",")).';
variants = split(strtrim(string(getenv("NSB_VARIANTS")))).';
if variants == "", variants = ["euler", "bdf2", "bdf2_exact_start"]; end
ladders = split(strtrim(string(getenv("NSB_LADDERS")))).';
modes = split(strtrim(string(getenv("NSB_MODES")))).';
if modes == "", modes = ["iterated", "semi-implicit"]; end
if ladders == "", ladders = ["flat", "circle"]; end
names = ["u_F", "sigma_F", "eta", "v", "sigma_P", "q", "p"];
primal = [true false true false false false true];
base = {"storage", 1.0, "friction_bjs", 1.0, "density_fluid", 1, "density_biot", 1, ...
    "convection", true};
n_fail = 0;

%% (a) bookkeeping
if contains(parts, "a")
    k = 2; h = 1/4; dt = 0.1; t0 = 0.5;
    f = patch_fields(k, "linear");
    mms = StokesBiotCoupling.get_manufactured_solution(["0"; "1"], f{:}, base{:});
    G = geometry("flat", h, 0);
    G.fm = set_tau(G.fm, mms, k, [t0, t0 + 3*dt]);
    for scheme = ["euler", "bdf2"]
        [fp, bp] = make_problems(k, mms, dt, scheme, t0);
        s0 = StokesBiotCoupling.initial_state(fp, bp, G.fm, G.bm, Order = k);
        s0.fluid.u_star = project_star(G.fm, mms.at(t0).u_F, k);
        d = mms.at(t0 + dt);
        [fp, bp] = set_data(fp, bp, d);
        [fs, bs, s1, info] = StokesBiotCoupling.solve_navier_stokes(fp, bp, G.fm, G.bm, G.curve, d, ...
            Order = k, friction_bjs = 1, state = s0, mode = "semi-implicit");
        [fs_m, bs_m] = manual_step(fp, bp, G, d, k, s0, s0.fluid.u_star);
        ok = isequal(fs.Unknowns.field(2).dof, fs_m.Unknowns.field(2).dof) ...
            && isequal(fs.Unknowns.field(1).dof, fs_m.Unknowns.field(1).dof) ...
            && isequal(bs.Unknowns.field(2).dof, bs_m.Unknowns.field(2).dof);
        n_fail = n_fail + report(ok, scheme + ": semi-implicit step = manual interface_problem + solve with beta = u*^n");
        u_star = divergence_free_velocity(fs);
        ok = abs(s1.time - (t0 + dt)) < 1e-15 && s1.step == 1 && info.iterations == 1 ...
            && ~info.stokes_presolve && isequal(s1.fluid.u_star, u_star) ...
            && isequal(info.u_star, u_star) && isequal(s1.fluid.u_star_previous, s0.fluid.u_star) ...
            && isequal(s1.fluid.velocity, fs.Unknowns.field(2).dof) ...
            && isequal(s1.fluid.velocity_previous, s0.fluid.velocity) ...
            && isequal(s1.biot.displacement, bs.Unknowns.field(2).dof);
        n_fail = n_fail + report(ok, scheme + ": state advanced once (time, step, u_star <- u* of fs, u_star_previous <- u_star, fluid, Biot)");
        % second step: the predictor of the scheme from both levels
        d = mms.at(t0 + 2*dt);
        [fp, bp] = set_data(fp, bp, d);
        fs2 = StokesBiotCoupling.solve_navier_stokes(fp, bp, G.fm, G.bm, G.curve, d, ...
            Order = k, friction_bjs = 1, state = s1, mode = "semi-implicit");
        if scheme == "euler"
            beta = s1.fluid.u_star;
        else
            beta = {2*s1.fluid.u_star{1} - s1.fluid.u_star_previous{1}; ...
                    2*s1.fluid.u_star{2} - s1.fluid.u_star_previous{2}};
        end
        fs2_m = manual_step(fp, bp, G, d, k, s1, beta);
        ok = isequal(fs2.Unknowns.field(2).dof, fs2_m.Unknowns.field(2).dof);
        predictors = ["u*^n", "2u*^n - u*^{n-1}"];
        n_fail = n_fail + report(ok, scheme + ": second step predicts with " + ...
            predictors(1 + (scheme == "bdf2")));
    end

    % the first step without u_star: a Stokes pre-solve, u_star_previous = []
    [fp, bp] = make_problems(k, mms, dt, "bdf2", t0);
    s0 = StokesBiotCoupling.initial_state(fp, bp, G.fm, G.bm, Order = k);
    d = mms.at(t0 + dt);
    [fp, bp] = set_data(fp, bp, d);
    [fs, ~, s1, info] = StokesBiotCoupling.solve_navier_stokes(fp, bp, G.fm, G.bm, G.curve, d, ...
        Order = k, friction_bjs = 1, state = s0, mode = "semi-implicit");
    fs_stokes = manual_step(fp, bp, G, d, k, s0, {});
    fs_m = manual_step(fp, bp, G, d, k, s0, divergence_free_velocity(fs_stokes));
    ok = info.stokes_presolve && info.iterations == 1 && isfield(s1.fluid, "u_star_previous") ...
        && isempty(s1.fluid.u_star_previous) && isequal(fs.Unknowns.field(2).dof, fs_m.Unknowns.field(2).dof);
    n_fail = n_fail + report(ok, "first step without u_star: Stokes pre-solve, beta = its u*, u_star_previous = []");
    [~, h2] = fluid_step_history("bdf2", dt, s1.fluid);
    n_fail = n_fail + report(isequal(h2.convective_velocity, s1.fluid.u_star), ...
        "after the first step: the bdf2 predictor is u*^1 (u_star_previous = [])");

    % Cache: iterated with and without the cache, then a step with info.Cache
    [fs_a, bs_a, s_a, info_a] = StokesBiotCoupling.solve_navier_stokes(fp, bp, G.fm, G.bm, G.curve, d, ...
        Order = k, friction_bjs = 1, state = s0, tolerance = 1e-10);
    [~, ~, ~, cache] = StokesBiotCoupling.interface_problem(fp, bp, G.fm, G.bm, G.curve, d, Order = k, ...
        friction_bjs = 1, state = s0);
    [fs_b, bs_b, ~, info_b] = StokesBiotCoupling.solve_navier_stokes(fp, bp, G.fm, G.bm, G.curve, d, ...
        Order = k, friction_bjs = 1, state = s0, tolerance = 1e-10, Cache = cache);
    ok = isequal(fs_a.Unknowns.field(2).dof, fs_b.Unknowns.field(2).dof) ...
        && isequal(bs_a.Unknowns.field(2).dof, bs_b.Unknowns.field(2).dof) ...
        && isequal(info_a.increments, info_b.increments) && ~isempty(info_a.Cache);
    n_fail = n_fail + report(ok, sprintf("Cache: identical iterated step (%d iterations, increments %s)", ...
        info_a.iterations, sprintf("%.1e ", info_a.increments)));
    d = mms.at(t0 + 2*dt);
    [fp, bp] = set_data(fp, bp, d);
    fs_c = StokesBiotCoupling.solve_navier_stokes(fp, bp, G.fm, G.bm, G.curve, d, ...
        Order = k, friction_bjs = 1, state = s_a, tolerance = 1e-10);
    fs_d = StokesBiotCoupling.solve_navier_stokes(fp, bp, G.fm, G.bm, G.curve, d, ...
        Order = k, friction_bjs = 1, state = s_a, tolerance = 1e-10, Cache = info_a.Cache);
    n_fail = n_fail + report(isequal(fs_c.Unknowns.field(2).dof, fs_d.Unknowns.field(2).dof), ...
        "Cache of the previous step (info.Cache): identical next step");
    % without convection, a Cache leaves interface_problem's solution unchanged
    [c1, ~, ~, cache0] = StokesBiotCoupling.interface_problem(fp, bp, G.fm, G.bm, G.curve, d, Order = k, ...
        friction_bjs = 1, state = s_a);
    c2 = StokesBiotCoupling.interface_problem(fp, bp, G.fm, G.bm, G.curve, d, Order = k, ...
        friction_bjs = 1, state = s_a, Cache = cache0);
    f1 = solve(c1, k);  f2 = solve(c2, k);
    n_fail = n_fail + report(isequal(f1.Unknowns.field(2).dof, f2.Unknowns.field(2).dof), ...
        "interface_problem without convection: Cache gives the identical solution");

    % the stabilization rule of interface_problem
    fq = fp;
    fq.Data.convective_velocity = s_a.fluid.u_star;
    fq.Data.convective_stabilization = true;
    n_fail = n_fail + report(throws(@() StokesBiotCoupling.interface_problem(fq, bp, G.fm, G.bm, ...
        G.curve, d, Order = k, friction_bjs = 1, state = s_a), "StokesBiotCoupling:stabilization"), ...
        "interface_problem: convective_velocity with convective_stabilization = true is an error");

    % s_margin, max_bn, max_bn_interface: the last solve's beta = u* = u_F
    % (polynomial of degree k < k + 1) at t0 + dt, evaluated independently
    ref = margin_reference(G.fm, mms.at(t0 + dt).u_F, k, 1, 1);
    dev = abs([info_a.s_margin, info_a.max_bn, info_a.max_bn_interface] - ...
        [ref.s_margin, ref.max_bn, ref.max_bn_interface]);
    n_fail = n_fail + report(all(dev < 1e-9), sprintf( ...
        "s_margin %.6f, max_bn %.6f, max_bn_interface %.6f = independent evaluation (dev %.1e)", ...
        info_a.s_margin, info_a.max_bn, info_a.max_bn_interface, max(dev)));
    fm_small = G.fm;
    fm_small.tau = 1e-3*ones(size(fm_small.tau));
    lastwarn("");
    warning("off", "StokesBiotCoupling:stabilization");
    [~, ~, ~, info] = StokesBiotCoupling.solve_navier_stokes(fp, bp, fm_small, G.bm, G.curve, ...
        mms.at(t0 + dt), Order = k, friction_bjs = 1, state = s0, mode = "semi-implicit");
    [~, id] = lastwarn();
    warning("on", "StokesBiotCoupling:stabilization");
    n_fail = n_fail + report(info.s_margin <= 0 && id == "StokesBiotCoupling:stabilization", ...
        sprintf("tau = 1e-3: s_margin %.3f <= 0 and the warning StokesBiotCoupling:stabilization", info.s_margin));
end

%% (b) patch test
if contains(parts, "b")
    degrees_b = degrees;
    if all(isnan(degrees_b)), degrees_b = 1:3; end
    tol = 1e-10;
    cases = struct("problem", {"flat", "flat", "circle"}, "gap", {"zero", "h2", "h2"}, ...
        "hs", {[1/4, 1/6], [1/4, 1/6], [1/8, 1/12]}, "asserted", {true, true, false});
    fprintf("\n(b) patch test, euler, linear in t, 3 steps from t = 1/2, dt = 0.1, iterated to 1e-12\n");
    for k = degrees_b
        f = patch_fields(k, "linear");
        mms = StokesBiotCoupling.get_manufactured_solution(["0"; "1"], f{:}, base{:});
        for c = cases
            if c.problem == "circle"
                mms_c = StokesBiotCoupling.get_manufactured_solution(circle_normal(), ...
                    f{:}, base{:});
            else
                mms_c = mms;
            end
            worst = 0;
            for h = c.hs
                G = geometry(c.problem, h, gap_of(c.gap, h));
                G.fm = set_tau(G.fm, mms_c, k, [0.5, 0.8]);
                out = march(k, G, mms_c, "euler", 0.1, 0.5, 3, "initial", "iterated", 1e-12);
                [w, at] = max(out.relative(end, :));
                worst = max(worst, w);
                fprintf("  k %d %-6s gap %-4s h 1/%-2d max rel. error %.2e (%s), Picard iterations %s, s_margin %.2f\n", ...
                    k, c.problem, c.gap, round(1/h), w, names(at), mat2str(out.iterations), min(out.s_margin));
                n_fail = n_fail + report(min(out.s_margin) > 0, "s_margin > 0", true);
            end
            if c.asserted
                n_fail = n_fail + report(worst <= tol, sprintf("(b) patch k = %d, %s gap %s: %.2e", ...
                    k, c.problem, c.gap, worst));
            else
                fprintf("  (reported) circle patch k = %d: %.2e\n", k, worst);
            end
        end
    end
end

%% (c) V5
if contains(parts, "c")
    k = 2; h = 1/8; dt = 1/20; t0 = 0.5;
    f = quadratic_smooth_fields();
    mms = StokesBiotCoupling.get_manufactured_solution(["0"; "1"], f{:}, base{:});
    G = geometry("flat", h, h^2);
    G.fm = set_tau(G.fm, mms, k, [t0, t0 + 4*dt]);
    out = march(k, G, mms, "bdf2", dt, t0, 3, "exact", "iterated", 1e-10);
    ratio = out.absolute(end, :)./out.absolute(1, :);
    fprintf("\n(c) V5: bdf2, exact start, k = %d, h = 1/%d, gap h^2, dt = %g, iterated (iterations %s)\n", ...
        k, round(1/h), dt, mat2str(out.iterations));
    print_table("step", 1:3, out.absolute, [], names);
    fprintf("  error(step 3)/error(step 1): %s\n", strjoin(names + " " + compose("%.2f", ratio), ", "));
    n_fail = n_fail + report(all(ratio <= 2) && min(out.s_margin) > 0, ...
        "(c) V5: every field's error after 3 bdf2 steps within 2x the first step's");
end

%% (d) V6 time order and (e) Picard vs semi-implicit
if contains(parts, "d")
    k = 2; h = 1/4; T = 1;
    dts = 1./[4, 8, 16, 32, 64];
    thresholds = struct("euler", 0.9, "bdf2", -Inf, "bdf2_exact_start", 1.8);
    mms = StokesBiotCoupling.get_manufactured_solution(["0"; "1"], base{:}, ...
        stream_function = "(2 + cos(2*t))*(x^3 - 2*x*y^2 + x^2*y + y^3/3 + x*y)", ...
        pressure_fluid = "(1 + sin(t))*(x^2 - x*y + 2*y^2 + x)", ...
        displacement = ["exp(t)*(x^2 + x*y + y)"; "(2 - cos(3*t))*(x*y - y^2 + x)"], ...
        pressure_biot = "cos(2*t)*(x^2 - x*y + y^2 + 1)");
    G = geometry("flat", h, h^2);
    G.fm = set_tau(G.fm, mms, k, linspace(0, T, 9));
    finals = struct();
    for variant = variants
        scheme = extractBefore(variant + "_", "_");
        start = "initial";
        if variant == "bdf2_exact_start", start = "exact"; end
        for mode = modes
            errors = zeros(numel(dts), numel(names));
            iters = zeros(numel(dts), 2);
            margin = Inf;
            last = cell(1, numel(dts));
            for i = 1:numel(dts)
                n_steps = round(T/dts(i));
                if start == "exact", n_steps = n_steps - 1; end
                out = march(k, G, mms, scheme, dts(i), 0, n_steps, start, mode, 1e-10);
                errors(i, :) = out.absolute(end, :);
                iters(i, :) = [mean(out.iterations), max(out.iterations)];
                margin = min(margin, min(out.s_margin));
                last{i} = out.fields;
            end
            finals.(variant).(strrep(mode, "-", "_")) = last;
            eoc = log2(errors(1:end-1, :)./errors(2:end, :));
            fprintf("\n(d) %s, %s: errors at T = %g, k = %d, h = 1/%d, gap h^2 (s_margin >= %.2f)\n", ...
                strrep(variant, "_", " "), mode, T, k, round(1/h), margin);
            print_table("dt", dts, errors, eoc, names);
            if mode == "iterated"
                fprintf("  Picard iterations per step (mean, max) per dt: %s\n", ...
                    strjoin(compose("%.2f/%d", iters(:, 1), iters(:, 2)), ", "));
            end
            threshold = thresholds.(variant);
            label = sprintf("(d) %s, %s: EOC in time of all fields (two finest steps) >= %.1f", ...
                strrep(variant, "_", " "), mode, threshold);
            if isfinite(threshold)
                n_fail = n_fail + report(all(eoc(end, :) >= threshold) && margin > 0, label);
            else
                fprintf("  (reported) %s, %s: EOC of the two finest steps %s\n", strrep(variant, "_", " "), ...
                    mode, strjoin(names + " " + compose("%.2f", eoc(end, :)), ", "));
            end
        end
    end
    if contains(parts, "e") && numel(modes) == 2
        for variant = variants
            it = finals.(variant).iterated;  si = finals.(variant).semi_implicit;
            D = zeros(numel(dts), numel(names));
            for i = 1:numel(dts)
                D(i, :) = field_differences(it{i}, si{i});
            end
            fprintf("\n(e) %s: || iterated - semi-implicit || at T (L2), rates\n", strrep(variant, "_", " "));
            print_table("dt", dts, D, log2(D(1:end-1, :)./D(2:end, :)), names);
        end
    end
end

%% (f) V7 space
if contains(parts, "f")
    degrees_f = degrees;
    if all(isnan(degrees_f)), degrees_f = 1:2; end
    dt = 0.2; t0 = 0.5;
    specs = struct("flat", struct("hs", [1/4, 1/8, 1/16], "slack", [0.3, 0.5]), ...
                   "circle", struct("hs", [1/10, 1/20, 1/40], "slack", [0.4, 0.6]));
    for problem = ladders
        S = specs.(problem);
        normal = ["0"; "1"];
        if problem == "circle", normal = circle_normal(); end
        f = quadratic_smooth_fields();
        mms = StokesBiotCoupling.get_manufactured_solution(normal, f{:}, base{:});
        for k = degrees_f
            errors = zeros(numel(S.hs), numel(names));
            margin = Inf;
            for i = 1:numel(S.hs)
                G = geometry(problem, S.hs(i), S.hs(i)^2);
                G.fm = set_tau(G.fm, mms, k, [t0, t0 + 3*dt]);
                out = march(k, G, mms, "bdf2", dt, t0, 2, "exact", "iterated", 1e-10);
                errors(i, :) = out.absolute(end, :);
                margin = min(margin, min(out.s_margin));
                fprintf("    %-6s k %d h 1/%-3d iterations %s s_margin %.2f max|div u*| %.1e\n", problem, k, ...
                    round(1/S.hs(i)), mat2str(out.iterations), min(out.s_margin), max(out.divergence));
            end
            rates = log(errors(1:end-1, :)./errors(2:end, :))./log(S.hs(1:end-1)./S.hs(2:end)).';
            fprintf("\n(f) V7 %s, k = %d, gap h^2, 2 bdf2 steps (exact start, quadratic in t), dt = %g\n", ...
                problem, k, dt);
            print_table("1/h", round(1./S.hs), errors, rates, names);
            overall = log(errors(1, :)./errors(end, :))/log(S.hs(1)/S.hs(end));
            need = k + 1 - S.slack(2)*ones(1, numel(names));
            need(primal) = k + 1 - S.slack(1);
            for j = 1:numel(names)
                n_fail = n_fail + report(overall(j) >= need(j), sprintf("(f) %s k = %d: %s rate %.2f >= %.1f", ...
                    problem, k, names(j), overall(j), need(j)));
            end
            n_fail = n_fail + report(margin > 0, sprintf("(f) %s k = %d: s_margin > 0", problem, k), true);
        end
    end
end

fprintf("\n%d check(s) failed.\n", n_fail);

%% local functions
function n = report(ok, label, quiet)
% print ok/FAIL and count the failure; quiet prints only failures
if nargin < 3, quiet = false; end
n = double(~ok);
if ~ok
    fprintf("  FAIL  %s\n", label);
elseif ~quiet
    fprintf("  ok    %s\n", label);
end
end

function f = patch_fields(k, time)
% polynomial fields in space: stream function of degree min(k, 2) + 1
% (u_F of degree <= 2 keeps (u . grad) u . v within the source rule at k
% = 3), the others of degree k, times a factor linear in t
g = ["(1 + 2*t)", "(2 - t)", "(1 + t)", "(3 - 2*t)", "(1 + 3*t)"];
if time == "quadratic"
    g = ["(1 + 2*t + t^2)", "(2 - t + 2*t^2)", "(1 + t - t^2)", "(3 - 2*t + t^2)", "(1 + 3*t + 2*t^2)"];
end
f = {"stream_function", g(1) + "*(" + poly_string(min(k, 2) + 1, 1) + ")", ...
     "pressure_fluid",  g(2) + "*(" + poly_string(k, 2) + ")", ...
     "displacement",   [g(3) + "*(" + poly_string(k, 3) + ")"; g(4) + "*(" + poly_string(k, 4) + ")"], ...
     "pressure_biot",   g(5) + "*(" + poly_string(k, 5) + ")"};
end

function f = quadratic_smooth_fields()
% smooth in space, quadratic in t (bdf2 is exact in time for them)
f = {"stream_function", "(2 + t - t^2)*sin(pi*x)*sin(pi*y)", ...
     "pressure_fluid", "(1 + t^2)*sin(pi*x)*exp(y - 0.5)", ...
     "displacement", ["(1 + t + t^2/2)*x*sin(y)"; "(2 - t + 2*t^2)*y*cos(x)"], ...
     "pressure_biot", "(1 - 2*t^2)*cos(pi*x)*sin(pi*y)"};
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

function n = circle_normal()
r = "sqrt((x - 0.5)^2 + (y - 0.5)^2)";
n = ["(x - 0.5)/" + r; "(y - 0.5)/" + r];
end

function g = gap_of(kind, h)
g = 0;
if kind == "h2", g = h^2; end
if kind == "h", g = h; end
end

function G = geometry(problem, h, gap)
switch problem
    case "flat"
        G.curve = FlatCurveHorizontal(0.5);
        [G.fm, G.bm] = fitted_split(G.curve, [0.8*h, h], gap);
    case "circle"
        G.curve = CircularCurve([0.5; 0.5], 0.25);
        [G.fm, G.bm] = interface_meshes(G.curve, "fitted", 0.8*h, gap = gap, h_solid = h);
end
end

function mesh = set_tau(mesh, mms, k, times)
% tau per face 1 + rho B_f / mu, B_f = max |u_F| over the face's Gauss
% points and the sample times (rho = mu = 1 here)
P = mms.parameters;
s = reshape(lgwt(k + 3, 0, 1), [], 1);
F = mesh.boundaries;
A = mesh.coordinates(1:2, mesh.faces(1, F(:)));
B = mesh.coordinates(1:2, mesh.faces(2, F(:)));
X = A(1, :) + s.*(B(1, :) - A(1, :));
Y = A(2, :) + s.*(B(2, :) - A(2, :));
bound = zeros(1, numel(F));
for t = linspace(min(times), max(times), 5)
    u = mms.at(t).u_F;
    bound = max(bound, max(sqrt(u{1}(X, Y).^2 + u{2}(X, Y).^2), [], 1));
end
n_faces = size(mesh.faces, 2);
face_max = accumarray(F(:), bound(:), [n_faces, 1], @max).';
mesh.tau = 1 + P.density_fluid*face_max(F)/P.viscosity;
end

function u_star = project_star(mesh, u, k)
% {u1; u2} projected onto P_{k+1} (the form of u*)
p1 = project(Datum(u{1}), mesh, DubinerBasis(k + 1));
p2 = project(Datum(u{2}), mesh, DubinerBasis(k + 1));
u_star = {p1.dof; p2.dof};
end

function ref = margin_reference(mesh, u, k, mu, rho)
% s_margin, max_bn, max_bn_interface of the exact u at the Gauss points
% of the element's convective face term (kq + 1 points, kq =
% ceil((k + 1 + 2k)/2)), with the exterior unit normals
kq = max(ceil((k + 1 + 2*k)/2), 1);
s = reshape(lgwt(kq + 1, 0, 1), [], 1);
F = mesh.boundaries;
A = mesh.coordinates(1:2, mesh.faces(1, F(:)));
B = mesh.coordinates(1:2, mesh.faces(2, F(:)));
N = reshape(mesh.exterior_unit_normals, 2, []);
X = A(1, :) + s.*(B(1, :) - A(1, :));
Y = A(2, :) + s.*(B(2, :) - A(2, :));
bn = N(1, :).*u{1}(X, Y) + N(2, :).*u{2}(X, Y);
ref.max_bn = max(abs(bn), [], "all");
ref.max_bn_interface = max(abs(bn(:, ismember(F(:), mesh.intfaces))), [], "all");
ref.s_margin = min(mu*mesh.tau(:).' - rho*max(bn, [], 1)/2);
end

function [fp, bp] = make_problems(k, mms, dt, scheme, t)
% the fluid and Biot problems with the initial data at time t
P = mms.parameters;
d = mms.at(t);
fp = IncompressibleElasticityProblem(k, viscosity = P.viscosity, ...
    density = P.density_fluid, time_step = dt, source_term = d.f_F, ...
    datum_dirichlet = d.u_F, pressure_mean = "pinned", initial_velocity = d.u_F, ...
    convective_stabilization = false);
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

function [fs, bs] = manual_step(fp, bp, G, d, k, state, beta)
% the coupled step with a given beta, by hand
fp.Data.convective_velocity = beta;
fp.Data.convective_stabilization = false;
coupled = StokesBiotCoupling.interface_problem(fp, bp, G.fm, G.bm, G.curve, d, Order = k, ...
    friction_bjs = 1, state = state);
[fs, bs] = solve(coupled, k);
end

function out = march(k, G, mms, scheme, dt, t0, n_steps, start, mode, tolerance)
% n_steps Navier-Stokes-Biot steps from t0 (start = "initial": the
% projected state at t0, no u_star (Stokes pre-solve at the first step),
% a bdf2 march starting with its Euler step; "exact": the projected
% states at t0 and t0 + dt as levels n-1 and n, with u_star,
% u_star_previous the P_{k+1} projections of u_F there, the march from t0
% + dt). out.absolute / out.relative: errors after each step (rows)
[fp, bp] = make_problems(k, mms, dt, scheme, t0);
state = StokesBiotCoupling.initial_state(fp, bp, G.fm, G.bm, Order = k);
if start == "exact"
    previous = state;
    [fp, bp] = make_problems(k, mms, dt, scheme, t0 + dt);
    state = StokesBiotCoupling.initial_state(fp, bp, G.fm, G.bm, Order = k);
    for name = ["displacement_skeleton", "stress", "rotation", "pressure", ...
            "displacement", "velocity"]
        state.biot.(name + "_previous") = previous.biot.(name);
    end
    state.fluid.velocity_previous = previous.fluid.velocity;
    state.fluid.u_star = project_star(G.fm, mms.at(t0 + dt).u_F, k);
    state.fluid.u_star_previous = project_star(G.fm, mms.at(t0).u_F, k);
    state.step = 1;
end
[out.absolute, out.relative] = deal(zeros(n_steps, 7));
[out.iterations, out.s_margin, out.divergence] = deal(zeros(1, n_steps));
cache = struct.empty;
for n = 1:n_steps
    d = mms.at(state.time + dt);
    [fp, bp] = set_data(fp, bp, d);
    [fs, bs, state, info] = StokesBiotCoupling.solve_navier_stokes(fp, bp, G.fm, G.bm, G.curve, d, ...
        Order = k, friction_bjs = mms.parameters.friction_bjs, state = state, mode = mode, ...
        tolerance = tolerance, max_iterations = 30, Cache = cache);
    cache = info.Cache;
    out.iterations(n) = info.iterations;
    out.s_margin(n) = info.s_margin;
    out.divergence(n) = info.divergence;
    [out.absolute(n, :), out.relative(n, :), out.fields] = field_errors(fs, bs, state, G.fm, G.bm, ...
        mms.at(state.time));
end
end

function [e, rel, list] = field_errors(fs, bs, state, fm, bm, d)
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

function D = field_differences(a, b)
% L2 norms of the differences of the fields of two runs (lists of
% field_errors)
D = zeros(1, size(a, 1));
for i = 1:size(a, 1)
    z = a{i, 1};
    z.dof = a{i, 1}.dof - b{i, 1}.dof;
    zero = @(x, y) 0*x;
    if iscell(a{i, 3}), zero = repmat({zero}, size(a{i, 3})); end
    D(i) = compute_l2_error(z, a{i, 2}, Datum(zero));
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
