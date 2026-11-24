% test_oseen_stress_ladder_l1.m
% Verification ladder, rung L1, for the stress-form Oseen fluid (PEERS
% stress sigma, velocity u, rotation gamma, with a prescribed convective
% velocity beta; the fluid of the Navier-Stokes-Biot coupling) alone on the
% unit square, k = 1, 2, 3 (k = 0 is skipped for the stress elements), on
% the configurations
%   matrix-pinned:      FEProblem + IncompressibleElasticityMatrices
%                       (convective_velocity = beta,
%                       convective_stabilization = false) + AffineOperator
%                       + ComponentwiseBoundaryConditions, traction face
%                       x = 1 (datum sigma n, n = (1, 0)), pressure_mean
%                       "pinned": the fluid as the coupling builds it;
%   matrix-constrained: the same, all Dirichlet, zero-mean pressure;
%   class-pinned:       solve(NavierStokesProblem, max_iterations = 1)
%                       with the same beta, convective_stabilization =
%                       false and the same mesh tau (a cross-check of the
%                       matrix path; one step, since the solve is a Picard
%                       loop whose later iterates take beta = u*).
% Unless stated otherwise the stabilisation is on the mesh, S = mu*tau with
% the face-wise rule tau = 1 + rho_c max_F |beta.n|/(2 mu) (rule "S"), so
% that S - rho_c beta.n/2 > 0 on every side of every element.
% beta cases: "u" = the exact velocity (handles; the manufactured
% Navier-Stokes field), "B" = (1 + y, x), "C" = -(1 + y, x); all are
% divergence free, so the conservative-form convective term of the scheme
% equals rho_c (beta.grad) u of the manufactured source. On the traction
% face x = 1, beta.n changes sign for "u" (inflow and outflow through the
% do-nothing face), is > 0 for "B" (outflow) and < 0 for "C" (inflow).
%  1. Patch test: stream function of degree k + 1 (u, p of degree <= k),
%     beta = "u", "B" as handles and as Dubiner dofs of degree k + 1 (the
%     form of u*: the L2 projection of the exact beta), steady, plus beta =
%     "u" transient (density 1, theta = 1, three steps, u linear in t), on
%     three meshes: every field reproduced to a relative L2 error <= 1e-10.
%     Exception (documented, reported, not counted): beta = "u" at k = 3.
%     There f = (beta.grad) u - div sigma has degree 2k - 1 = 5 and
%     test_elements(mesh, f, k) integrates f v (degree 3k - 1 = 8) with the
%     rule of order k (exact to degree 2k + 1 = 7), so the load vector is
%     not exact and the patch error is the quadrature error of the source.
%  2. Steady convergence, smooth solution, n = 4..32: last EOC >= k + 1 -
%     0.3 for sigma, u, gamma, p, mu = 0.1, beta = "u", "B" (all three
%     configurations) and "C" (matrix-pinned). Convection-dominated case
%     mu = 1e-3, beta = "u", "B" (matrix-pinned): reported only.
%  3. Transient Oseen, density 1, theta = 1, three steps of dt = 0.1 from
%     Pi u(0), beta = u(t_{n+1}) (handles), u linear in t (time error zero):
%     same rates, all three configurations.
%  4. Stabilisation, matrix-pinned, beta = "u", mu = 0.1:
%     (a) convective_stabilization = true, mesh tau = 1 (the class rule
%         S = mu + rho_c max|beta.n|/2): asserted;
%     (b) false with rule "S" (as in 2): asserted;
%     (c) false with mesh tau = 1 and tau = 1e-2, which violate
%         S - rho_c beta.n/2 > 0 on some sides: reported only (min margin,
%         share of violating sides, errors, rates);
%     (d) mesh tau = 500 (S = 500 mu, the production level of S; it
%         satisfies the S-rule here), k = 2, error constants relative to
%         (b) reported. Documented exception, reported only: the gamma
%         rate is about 2.3 at n = 32 and rises only slowly (2.51, 2.67 at
%         n = 64, 128), sigma about 2.75, and the sigma, gamma, p
%         constants are 30-50 times those of (b). The same numbers come
%         out with beta = 0 (plain Stokes, tau = 500 against tau = 1): a
%         property of the stress-form element at tau h >> 1, not of the
%         convective term;
%     (d') rule "S" times 500 (tau up to about 8000, S/mu far above the
%         production level), k = 2: reported only. Its gamma rate drops to
%         about 2.2-2.3 and its sigma, p constants grow by 150-290 over (b)
%         on n = 4..128 (pre-asymptotic: tau h >> 1 on every level).
% Pieces: set OSEEN_L1_PARTS to a comma list of patch, steady, transient,
% stab (default all) and OSEEN_L1_K to a comma list of degrees (default
% 1,2,3), e.g. OSEEN_L1_PARTS=patch OSEEN_L1_K=2; OSEEN_L1_NS overrides the
% levels of 2-4 (default 4,8,16,32).
clc; clear;
n_fail = 0;

parts = string(getenv("OSEEN_L1_PARTS"));
if parts == "", parts = "patch,steady,transient,stab"; end
parts = strtrim(split(parts, ",")).';
ks = str2double(split(string(getenv("OSEEN_L1_K")), ",")).';
if all(isnan(ks)), ks = 1:3; end

mu = 0.1; rho = 1; dt = 0.1; n_steps = 3;
configs = ["matrix-pinned", "matrix-constrained", "class-pinned"];
ns = str2double(split(string(getenv("OSEEN_L1_NS")), ",")).';
if all(isnan(ns)), ns = [4 8 16 32]; end
U = "(1 + t)*sin(pi*x + 0.3)*sin(pi*y + 0.7)";
P = "(1 + t^2)*cos(pi*x)*sin(pi*y) + x*y";
betas = struct("u", [], "B", ["1 + y", "x"], "C", ["-1 - y", "-x"]);

%% 1. patch test
if ismember("patch", parts)
    patch_U = ["x^2 + 3*x*y - y^2 + 2*x", ...
               "x^3 + 2*x^2*y - x*y^2 + y^3 + x*y", ...
               "x^4 - 2*x^3*y + x^2*y^2 + y^4 + x^2*y"];
    patch_P = ["x - 2*y", "x^2 - x*y + 3*y", "x^3 - x*y^2 + y^2 + x"];
    patch_meshes = {[2 2], [3 5], [6 4]};
    for k = ks
        Up = "(1 + t)*(" + patch_U(k) + ")";
        Pp = "(1 + t^2)*(" + patch_P(k) + ")";
        cases = ["u", "u-dofs", "B", "B-dofs", "u-transient"];
        for c = cases
            bname = extractBefore(c + "-", "-");
            transient = endsWith(c, "transient");
            if transient
                d_at = arrayfun(@(n) NavierStokesProblem.get_manufactured_solution( ...
                    Up, Pp, mu, betas.(bname), rho, n*dt), 0:n_steps);
            else
                d_at = NavierStokesProblem.get_manufactured_solution(Up, Pp, mu, betas.(bname), 0, 0);
            end
            documented = bname == "u" && k == 3;
            e_max = 0;
            for cfg = configs
                for m = 1:numel(patch_meshes)
                    nm = patch_meshes{m};
                    mesh = StructuredSimplicialGrid(nm(1), nm(2));
                    [sol, d] = run_case(cfg, k, mesh, d_at, mu, rho, dt, ...
                        endsWith(c, "dofs"), false, 1);
                    e = field_errors(sol, mesh, k, d) ./ field_scales(mesh, d);
                    e_max = max(e_max, max(e));
                    fprintf("patch %-12s %-18s k = %d mesh %dx%d: rel. err sigma %.1e  u %.1e  gamma %.1e  p %.1e\n", ...
                        c, cfg, k, nm, e);
                    if max(e) > 1e-10 && ~documented
                        fprintf("patch %s %s k = %d mesh %dx%d FAILED\n", c, cfg, k, nm);
                        n_fail = n_fail + 1;
                    end
                end
            end
            fprintf("patch %-12s k = %d: max relative error %.2e%s\n", c, k, e_max, ...
                ternary(documented, "  (documented: source quadrature, see header)", ""));
        end
    end
end

%% 2. steady convergence
if ismember("steady", parts)
    for k = ks
        for bname = ["u", "B", "C"]
            d = NavierStokesProblem.get_manufactured_solution(U, P, mu, betas.(bname), 0, 0);
            cfgs = configs;
            if bname == "C", cfgs = "matrix-pinned"; end
            for cfg = cfgs
                E = rates_study(cfg, k, ns, d, mu, rho, dt, "S", 1);
                n_fail = n_fail + report(sprintf("steady mu = 0.1 beta = %s %s k = %d", bname, cfg, k), ns, E, k + 1, true);
            end
        end
        for bname = ["u", "B"]
            d = NavierStokesProblem.get_manufactured_solution(U, P, 1e-3, betas.(bname), 0, 0);
            E = rates_study("matrix-pinned", k, ns, d, 1e-3, rho, dt, "S", 1);
            report(sprintf("steady mu = 1e-3 beta = %s matrix-pinned k = %d", bname, k), ns, E, k + 1, false);
        end
    end
end

%% 3. transient convergence
if ismember("transient", parts)
    d_at = arrayfun(@(n) NavierStokesProblem.get_manufactured_solution( ...
        U, P, mu, [], rho, n*dt), 0:n_steps);
    for k = ks
        for cfg = configs
            E = rates_study(cfg, k, ns, d_at, mu, rho, dt, "S", 1);
            n_fail = n_fail + report(sprintf("transient beta = u %s k = %d", cfg, k), ns, E, k + 1, true);
        end
    end
end

%% 4. stabilisation
if ismember("stab", parts)
    d = NavierStokesProblem.get_manufactured_solution(U, P, mu, [], 0, 0);
    for k = ks
        E = rates_study("matrix-pinned", k, ns, d, mu, rho, dt, "class", 1);
        n_fail = n_fail + report(sprintf("stab (a) class rule k = %d", k), ns, E, k + 1, true);
        Eb = rates_study("matrix-pinned", k, ns, d, mu, rho, dt, "S", 1);
        n_fail = n_fail + report(sprintf("stab (b) rule S k = %d", k), ns, Eb, k + 1, true);
        fprintf("stab (a)/(b) error ratios k = %d (sigma u gamma p, per level):\n", k);
        fprintf("  %3d   %.3f  %.3f  %.3f  %.3f\n", [ns; (E ./ Eb).']);
        for tau_c = [1, 1e-2]
            [E, margin] = rates_study("matrix-pinned", k, ns, d, mu, rho, dt, "fixed", tau_c);
            report(sprintf("stab (c) tau = %g k = %d", tau_c, k), ns, E, k + 1, false);
            fprintf("stab (c) tau = %g k = %d: min (S - rho beta.n/2)/mu per level %s; violating sides %s\n", ...
                tau_c, k, join(compose("%.3g", margin(:, 1).'), " / "), ...
                join(compose("%.1f%%", 100*margin(:, 2).'), " / "));
        end
        if k == 2
            Ed = rates_study("matrix-pinned", k, ns, d, mu, rho, dt, "fixed", 500);
            report("stab (d) tau = 500 k = 2", ns, Ed, k + 1, false);
            fprintf("stab (d)/(b) error ratios k = 2 (sigma u gamma p, per level):\n");
            fprintf("  %3d   %.3f  %.3f  %.3f  %.3f\n", [ns; (Ed ./ Eb).']);
            Ed = rates_study("matrix-pinned", k, ns, d, mu, rho, dt, "S", 500);
            report("stab (d') rule S x 500 k = 2", ns, Ed, k + 1, false);
            fprintf("stab (d')/(b) error ratios k = 2 (sigma u gamma p, per level):\n");
            fprintf("  %3d   %.3f  %.3f  %.3f  %.3f\n", [ns; (Ed ./ Eb).']);
        end
    end
end

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_oseen_stress_ladder_l1:failed", "%d checks failed", n_fail);
end

%% ------------------------------------------------------------------------
function [E, margin] = rates_study(cfg, k, ns, d_at, mu, rho, dt, rule, factor)
% errors [sigma, u, gamma, p] per level and, per level, [min (S - rho_c
% beta.n/2)/mu over sides and face points, share of sides where it is <= 0]
E = zeros(numel(ns), 4);
margin = zeros(numel(ns), 2);
for i = 1:numel(ns)
    mesh = StructuredSimplicialGrid(ns(i), ns(i));
    [sol, d, margin(i, :)] = run_case(cfg, k, mesh, d_at, mu, rho, dt, false, rule, factor);
    E(i, :) = field_errors(sol, mesh, k, d);
end
end

function [sol, d_end, margin] = run_case(cfg, k, mesh, d_at, mu, rho, dt, as_dofs, rule, factor)
% one steady solve (numel(d_at) == 1) or backward Euler from u^0 = Pi u(0)
% (d_at(n + 1) holds the data at t = n*dt); beta = d_at(n + 1).beta, as
% handles or projected onto dubiner_basis(k + 1) (as_dofs). The mesh tau:
% rule "S" (times factor), "class" (tau = 1, convective_stabilization =
% true) or "fixed" (tau = factor); false (the patch test) is rule "S"
transient = numel(d_at) > 1;
rho_c = ternary(transient, rho, 1);
if isequal(rule, false), rule = "S"; end
conv_stab = rule == "class";
bn_max = side_beta_n(mesh, d_at(end).beta);   % tau fixed in time: from the last (largest) beta
switch rule
    case "S",     mesh.tau = factor*(1 + rho_c*max(abs(bn_max), [], 3)/(2*mu));
    case "class", mesh.tau = ones(size(mesh.tau));
    case "fixed", mesh.tau = factor*ones(size(mesh.tau));
end
bn_out = max(bn_max, [], 3);                  % max outward beta.n_K per side
if conv_stab
    m = (mu + rho_c*max(abs(bn_max), [], 'all')/2 - rho_c*bn_out/2)/mu;
else
    m = (mu*mesh.tau - rho_c*bn_out/2)/mu;
end
margin = [min(m, [], 'all'), mean(m(:) <= 0)];
d2 = nchoosek(k + 2, 2);
if ~transient
    sol = solve_config(cfg, k, mesh, d_at, mu, beta_of(mesh, d_at.beta, k, as_dofs), {}, struct.empty, conv_stab);
    d_end = d_at;
    return
end
u_prev = L2proj2D(mesh, d_at(1).u, k);
for n = 1:numel(d_at) - 1
    opts = {"density", rho, "time_step", dt, "theta", 1.0};
    state = struct("velocity", u_prev);
    sol = solve_config(cfg, k, mesh, d_at(n + 1), mu, beta_of(mesh, d_at(n + 1).beta, k, as_dofs), ...
        opts, state, conv_stab);
    u_prev = reshape(sol.Unknowns.field(2).dof, 2*d2, []);
end
d_end = d_at(end);
end

function b = beta_of(mesh, beta, k, as_dofs)
% beta as handles or as {bx; by} dofs in dubiner_basis(k + 1)
b = beta;
if as_dofs
    b = {L2proj2D(mesh, beta{1}, k + 1); L2proj2D(mesh, beta{2}, k + 1)};
end
end

function bn = side_beta_n(mesh, beta)
% beta.n_K (outward, unit normal) at 6 Gauss points of every side: 3 x
% n_elements x 6, the side numbering of mesh.tau and convective_matrix
n_el = size(mesh.elements, 2);
[s, ~] = get_1d_reference_quadrature(5);
nq = numel(s);
q2 = [[0*s(:), 1 - s(:)]; [1 - s(:), 0*s(:)]; [1 - s(:), s(:)]];
[F, xy_0] = get_affine_map(mesh, output = "array");
xy = pagemtimes(F, q2.') + reshape(xy_0, 2, 1, []);
x = reshape(xy(1, :, :), [], n_el);  y = reshape(xy(2, :, :), [], n_el);
bx = beta{1}(x, y) + 0*x;  by = beta{2}(x, y) + 0*x;
area = mesh.area(mesh.boundaries);
nrm = mesh.normals .* kron(mesh.orientation, [1; 1]);
bn = zeros(3, n_el, nq);
for e = 1:3
    rows = (e - 1)*nq + (1:nq);
    v = (nrm(2*e - 1, :) .* bx(rows, :) + nrm(2*e, :) .* by(rows, :)) ./ area(e, :);
    bn(e, :, :) = reshape(v.', 1, n_el, nq);
end
end

function sol = solve_config(cfg, k, mesh, d, mu, beta, opts, state, conv_stab)
% one solve of configuration cfg with the data d (source, Dirichlet datum
% u, traction datum sigma n on x = 1 under "pinned") and the convective
% velocity beta; opts and state are the transient options and the velocity
% history (empty: steady)
parts = split(cfg, "-");
[path, closure] = deal(parts(1), parts(2));
gN = {@(x, y) 0*x; @(x, y) 0*x};
if closure == "pinned"
    mesh = apply_tag(mesh, "Neumann", @(X) isapprox(X(1, :), 1));
    gN = {d.sigma{1,1}; d.sigma{2,1}};
end
if path == "matrix"
    Pk = DubinerBasis(k);
    gspaces = [LegendreBasis(k).to_vector(2), DubinerBasis(0)];
    if closure == "constrained"
        gspaces = [gspaces, ConstantBasis()];
    end
    pb = FEProblem(Order = k, Spaces = [ ...
        FESpace([PEERSBasis(k), Pk.to_vector(), Pk, DubinerBasis(0)], Condense = true), ...
        FESpace(gspaces)]);
    pb.Mesh = mesh;
    bc = ComponentwiseBoundaryConditions([ ...
        boundary_condition("Dirichlet", x_data = d.u{1}, y_data = d.u{2});
        boundary_condition("Neumann", x = "neumann", y = "neumann", ...
            x_data = gN{1}, y_data = gN{2})]);
    if ~isempty(state), opts = [opts, {"state", state}]; end
    [lhs, flux, rhs] = IncompressibleElasticityMatrices(k, mesh, mu, d.f, ...
        "pressure_mean", closure, "natural_faces", bc.natural_faces(mesh), ...
        "convective_velocity", beta, "convective_stabilization", conv_stab, opts{:});
    pb.Operator = AffineOperator(lhs + rhs) + flux;
    pb.BoundaryConditions = bc;
    sol = solve(pb, mesh, k);
else
    args = {"viscosity", mu, "pressure_mean", closure, "source_term", d.f, ...
        "datum_dirichlet", d.u, "datum_neumann", gN, "convective_velocity", beta, ...
        "convective_stabilization", conv_stab};
    if ~isempty(opts)
        args = [args, {"density", opts{2}, "time_step", opts{4}, "time_scheme", "euler", ...
            "state", struct("fluid", state)}];
    end
    % one Oseen solve: NavierStokesProblem.solve is a Picard loop whose
    % later iterates replace beta by u*
    sol = solve(NavierStokesProblem(k, args{:}), mesh, k, max_iterations = 1);
end
end

function E = field_errors(sol, mesh, k, d)
% L2 errors [sigma, u, gamma, p], p_h = -tr(sigma_h)/2
e = sol.compute_l2_error({d.sigma, d.u, d.gamma});
E = [e(1:3).', pressure_error(mesh, k, sol.Unknowns.field(1).dof, d.p)];
end

function e = pressure_error(mesh, k, sigma_dof, p_ex)
[xy, rs, w] = quadpoints_elements_2d(mesh, k + 2, output = "array");
B = PEERSBasis(k).evaluate(rs, mesh);
S = pagemtimes(B, reshape(sigma_dof, size(B, 2), 1, 1, 1, []));
p_h = -reshape(S(:, 1, 1, 1, :) + S(:, 1, 2, 2, :), size(S, 1), [])/2;
p = reshape(p_ex(xy(1, :, :), xy(2, :, :)) + 0*xy(1, :, :), size(p_h));
e = sqrt(sum(mesh.volume .* (w'*(p_h - p).^2)));
end

function s = field_scales(mesh, d)
% L2 norms [sigma, u, gamma, p] of the exact fields, the velocity norm
% standing in for a field that vanishes
s = [l2_norm(mesh, d.sigma), l2_norm(mesh, d.u), l2_norm(mesh, {d.gamma}), l2_norm(mesh, {d.p})];
s(s < 1e-12*s(2)) = s(2);
end

function n = l2_norm(mesh, fs)
[xy, ~, w] = quadpoints_elements_2d(mesh, 6, output = "array");
n = 0;
for i = 1:numel(fs)
    v = reshape(fs{i}(xy(1, :, :), xy(2, :, :)) + 0*xy(1, :, :), size(xy, 2), []);
    n = n + sum(mesh.volume .* (w'*v.^2));
end
n = sqrt(n);
end

function n_fail = report(name, ns, E, expected, asserted)
fprintf("%s\n   N    |sigma|     EOC    |u|         EOC    |gamma|     EOC    |p|         EOC\n", name);
eoc = compute_eoc(E, 1./ns);
for i = 1:numel(ns)
    r = NaN(1, 4);
    if i > 1, r = eoc(i-1, :); end
    fprintf("  %3d   %.3e  %5.2f   %.3e  %5.2f   %.3e  %5.2f   %.3e  %5.2f\n", ns(i), ...
        E(i, 1), r(1), E(i, 2), r(2), E(i, 3), r(3), E(i, 4), r(4));
end
n_fail = 0;
if asserted && any(~(r >= expected - 0.3))
    n_fail = 1;
    fprintf("%s: last EOC below k + 1 - 0.3 = %.1f FAILED\n", name, expected - 0.3);
elseif ~asserted
    fprintf("%s: reported only\n", name);
end
end

function v = ternary(c, a, b)
if c, v = a; else, v = b; end
end
