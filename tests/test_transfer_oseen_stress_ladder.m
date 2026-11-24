% test_transfer_oseen_stress_ladder.m
% Verification ladder, rung L2, for the Oseen fluid in stress form (PEERS
% stress sigma, velocity u, rotation gamma, p_h = -tr(sigma_h)/2): the
% equation alone on a non-conforming boundary, the Dirichlet velocity
% transferred along the paths (transfer_mapping ->
% apply_transfer_conditions, path gradient grad u = A sigma + gamma). The
% domain is the unit square; the mesh covers [0,1] x [0,1-delta] and the
% datum on y = 1 is transferred along vertical paths of length delta.
% mu = 0.1, k = 1, 2, 3 (k = 0 is skipped for the PEERS element).
%
% Path: the class path only, IncompressibleElasticityProblem with
% Data.convective_velocity set to the prescribed beta (one Oseen solve, as
% each step of NavierStokesProblem.solve; stage 1 checks that the two agree
% with transfer). The matrix path (IncompressibleElasticityMatrices) has no
% transfer rows for this fluid.
%
% Convective velocity, prescribed as exact handles:
%   "u":      beta = the exact velocity; on the transferred boundary y = 1
%             beta.n = u_y changes sign (inflow and outflow faces);
%   "affine": beta = (1 + y, x); beta.n = x >= 0 there (outflow only).
% Stabilization, two configurations:
%   "fixed":  convective_stabilization = false, S = mu*tau with tau =
%             1 + rho_c B/(2 mu) on the mesh, B = max|beta| over the unit
%             square, so S - rho_c beta.n/2 >= mu > 0 on every face (the
%             configuration of the coupling);
%   "class":  convective_stabilization = true, mesh tau = 1: the class rule
%             S = mu + rho_c max|beta.n|/2.
% On a transfer face beta.n is evaluated on the mesh face (Gamma_h), in
% the local solver's face term <(uhat x beta) n, v>; the face's global
% flux row, the only place of the convective flux BbLL, is replaced by
% the transfer row, so no convective flux enters there and nothing is
% evaluated at the transferred point.
% The exact pressure is shifted to zero mean over the strip (the discrete
% one has zero mean over the mesh); sigma's diagonal accordingly.
%
% Stage 1 (consistency): fitted square (delta = 0), identity mapping: equal
%   to the plain Dirichlet Oseen solve (max dof difference <= 1e-9), and
%   NavierStokesProblem.solve(max_iterations = 1) with transfer equals the
%   single Oseen solve; polynomial data, n = 4.
% Stage 2 (patch test): u and p polynomial of degree <= k (stream function
%   of degree k + 1), so sigma and gamma are in the discrete spaces;
%   delta = h^2 and delta = h, n = 4, 8, 16: relative L2 errors of
%   sigma, u, gamma, p <= max(1e-10, 100 r), r the largest round-off of
%   the three global solves measured by one step of iterative refinement
%   (solve_roundoff): with the gap h at k = 3 the error reaches 1.2e-10 at
%   n = 16 while staying at round-off. Not asserted for beta = u at k = 3: the
%   source term (beta.grad)u has degree 2k - 1, so f v has degree
%   3k - 1 = 8, and the element rule of the source (test_elements, order
%   k) is exact to degree 2k + 1 = 7 only: the patch is not reproduced, with or without transfer (fitted
%   square, plain Dirichlet, n = 4, 8, 16: e(sigma) 4.5e-6, 6.9e-8,
%   1.1e-9; on the strip, relative 4.2e-6, 7.3e-8, 1.2e-9 at delta = h^2).
% Stage 3 (rates): stream function sin(pi x + 0.3) sin(pi y + 0.7),
%   p = cos(pi x) sin(pi y) + x y, n = 4..32 (k = 1 also 64). delta = h^2:
%   the EOC of the last pair must be >= k + 1 - 0.3 for sigma, u, gamma
%   and p. delta = h: reported only.
% Stage 4 (transient): density 1, theta = 1 (backward Euler), dt = 0.1,
%   three steps from u^0 = Pi u(0), u linear in t (p quadratic), so Euler
%   is exact in time: the patch test of stage 2 (delta = h^2 and h) and the
%   rates of stage 3 (delta = h^2 asserted, delta = h reported) at t = 0.3.
%   The patch with delta = h at k = 3 is reported, not asserted: round-off
%   grows over the steps (beta affine: relative 8e-11, 1e-11, 5e-9 at
%   n = 4, 8, 16), the regime in which the k = 3 rates with delta = h
%   stagnate (steady, beta affine: e(sigma) 2.7e-4, 3.4e-4, 4.3e-4 at
%   n = 16, 32, 64; the Stokes fluid of test_transfer_bc_incompressible_
%   elasticity diverges there too).
%
% To run a part, set before calling: HDG_TEST_K (e.g. "2") and
% HDG_TEST_STAGES (e.g. "3 4"); both default to everything.
clc; clear;
n_fail = 0;
k_list = env_list("HDG_TEST_K", 1:3);
stages = env_list("HDG_TEST_STAGES", 1:4);

mu = 0.1;  rho = 1;  dt = 0.1;  n_steps = 3;
betas = ["u", "affine"];
configs = ["fixed", "class"];
names = ["sigma", "u", "gamma", "p"];
tol_patch = 1e-10;

U_smooth = "(1 + t)*sin(pi*x + 0.3)*sin(pi*y + 0.7)";
P_smooth = "(1 + t^2)*(cos(pi*x)*sin(pi*y) + x*y)";
U_patch = ["x^2 + 3*x*y - y^2 + 2*x", ...
           "x^3 + 2*x^2*y - x*y^2 + y^3 + x*y", ...
           "x^4 - 2*x^3*y + x^2*y^2 + y^4 + x^2*y"];
P_patch = ["x - 2*y", "x^2 - x*y + 3*y", "x^3 - x*y^2 + y^2 + x"];

%% Stage 1: consistency at delta = 0
if ismember(1, stages)
    mesh0 = StructuredSimplicialGrid(4, 4);
    for k = k_list
        for b = betas
            d = manufactured(patch_U(U_patch, k), patch_P(P_patch, k), mu, b, 0, 0);
            for cfg = configs
                mesh = with_tau(mesh0, cfg, d.beta, 1, mu);
                pb = oseen_problem(k, mu, d, cfg, 1);
                s0 = solve(pb, mesh, Order = k);
                pb.Data.transfer_mapping = @(X) X;
                s1 = solve(pb, mesh, Order = k);
                ns_pb = NavierStokesProblem(k, viscosity = mu, source_term = d.f, ...
                    datum_dirichlet = d.u, convective_velocity = d.beta, transfer_mapping = @(X) X, ...
                    convective_stabilization = cfg == "class");
                s2 = solve(ns_pb, mesh, k, max_iterations = 1);
                [e01, e12] = deal(0);
                for f = 1:3
                    e01 = max(e01, max(abs(s0.Unknowns.field(f).dof - s1.Unknowns.field(f).dof), [], 'all'));
                    e12 = max(e12, max(abs(s1.Unknowns.field(f).dof - s2.Unknowns.field(f).dof), [], 'all'));
                end
                fprintf("Stage 1 k = %d beta %-6s %-5s: identity transfer vs Dirichlet %.1e, " + ...
                    "NavierStokesProblem (1 it) vs Oseen solve %.1e\n", k, b, cfg, e01, e12);
                if ~(max(e01, e12) <= 1e-9)
                    fprintf("Stage 1 k = %d beta %s %s FAILED\n", k, b, cfg); n_fail = n_fail + 1;
                end
            end
        end
    end
end

%% Stage 2: patch test on the strip
if ismember(2, stages)
    for k = k_list
        for b = betas
            d = manufactured(patch_U(U_patch, k), patch_P(P_patch, k), mu, b, 0, 0);
            for cfg = configs
                for gap = ["h2", "h"]
                    [worst, roundoff] = deal(zeros(1, 3));
                    for i = 1:3
                        [sol, mesh, gs] = solve_strip(k, 2^(i + 1), gap, mu, d, cfg, 0, {}, struct.empty);
                        [~, r] = field_errors(sol, mesh, k, d);
                        worst(i) = max(r);
                        roundoff(i) = solve_roundoff(gs);
                    end
                    bound = max(tol_patch, 100*max(roundoff));
                    fprintf("Stage 2 k = %d beta %-6s %-5s gap %-2s: max rel err (sigma, u, gamma, p) n = 4, 8, 16: %.1e %.1e %.1e (bound %.1e)\n", ...
                        k, b, cfg, gap, worst, bound);
                    worst = max(worst);
                    if ~source_exact(k, b)
                        fprintf("Stage 2 k = %d beta %s %s gap %s: reported, not asserted (source quadrature, see header)\n", k, b, cfg, gap);
                    elseif ~(worst <= bound)
                        fprintf("Stage 2 k = %d beta %s %s gap %s FAILED\n", k, b, cfg, gap); n_fail = n_fail + 1;
                    end
                end
            end
        end
    end
end

%% Stage 3: smooth solution, rates
if ismember(3, stages)
    for k = k_list
        ns = [4 8 16 32];
        if k == 1, ns = [ns 64]; end %#ok<AGROW>
        for b = betas
            d = manufactured(U_smooth, P_smooth, mu, b, 0, 0);
            for cfg = configs
                for gap = ["h2", "h"]
                    E = zeros(numel(ns), 4);
                    for i = 1:numel(ns)
                        [sol, mesh] = solve_strip(k, ns(i), gap, mu, d, cfg, 0, {}, struct.empty);
                        E(i, :) = field_errors(sol, mesh, k, d);
                    end
                    name = sprintf("steady k = %d beta %s %s gap %s", k, b, cfg, gap);
                    n_fail = n_fail + report(name, ns, gap, E, k + 1 - 0.3, gap == "h2");
                end
            end
        end
    end
end

%% Stage 4: transient, backward Euler, three steps, linear in t
if ismember(4, stages)
    for k = k_list
        for b = betas
            dp = arrayfun(@(n) manufactured(patch_U(U_patch, k, true), patch_P(P_patch, k, true), ...
                mu, b, rho, n*dt), 0:n_steps);
            ds = arrayfun(@(n) manufactured(U_smooth, P_smooth, mu, b, rho, n*dt), 0:n_steps);
            for cfg = configs
                for gap = ["h2", "h"]
                    worst = zeros(1, 3);
                    for i = 1:3
                        [sol, mesh] = run_transient(k, 2^(i + 1), gap, mu, dp, cfg, rho, dt);
                        [~, r] = field_errors(sol, mesh, k, dp(end));
                        worst(i) = max(r);
                    end
                    fprintf("Stage 4 patch k = %d beta %-6s %-5s gap %-2s: max rel err at t = %.1f, n = 4, 8, 16: %.1e %.1e %.1e\n", ...
                        k, b, cfg, gap, n_steps*dt, worst);
                    worst = max(worst);
                    if ~source_exact(k, b)
                        fprintf("Stage 4 patch k = %d beta %s %s gap %s: reported, not asserted (source quadrature, see header)\n", k, b, cfg, gap);
                    elseif k == 3 && gap == "h"
                        fprintf("Stage 4 patch k = %d beta %s %s gap %s: reported, not asserted (O(h) gap at k = 3, see header)\n", k, b, cfg, gap);
                    elseif ~(worst <= tol_patch)
                        fprintf("Stage 4 patch k = %d beta %s %s gap %s FAILED\n", k, b, cfg, gap);
                        n_fail = n_fail + 1;
                    end
                end
                for gap = ["h2", "h"]
                    ns = [4 8 16 32];
                    E = zeros(numel(ns), 4);
                    for i = 1:numel(ns)
                        [sol, mesh] = run_transient(k, ns(i), gap, mu, ds, cfg, rho, dt);
                        E(i, :) = field_errors(sol, mesh, k, ds(end));
                    end
                    name = sprintf("transient k = %d beta %s %s gap %s", k, b, cfg, gap);
                    n_fail = n_fail + report(name, ns, gap, E, k + 1 - 0.3, gap == "h2");
                end
            end
        end
    end
end

fprintf("%d check(s) failed.\n", n_fail);

%% ------------------------------------------------------------------ %%
function v = env_list(name, default)
% numeric list from an environment variable ("1 2"), default when unset
v = default;
s = strtrim(getenv(name));
if ~isempty(s), v = str2num(s); end %#ok<ST2NM>
end

function tf = source_exact(k, b)
% the element rule of the source (test_elements, order k) is exact to
% degree 2k + 1; with beta = u of degree k the source f has degree 2k - 1,
% so f v is integrated exactly only for k <= 2
tf = b ~= "u" || 3*k - 1 <= 2*k + 1;
end

function U = patch_U(U_patch, k, transient)
U = U_patch(k);
if nargin > 2 && transient, U = "(1 + t)*(" + U + ")"; end
end

function P = patch_P(P_patch, k, transient)
P = P_patch(k);
if nargin > 2 && transient, P = "(1 + t^2)*(" + P + ")"; end
end

function d = manufactured(U, P, mu, b, rho, t)
% beta "u": the exact velocity; "affine": (1 + y, x)
if b == "u", beta = []; else, beta = ["1 + y", "x"]; end
d = NavierStokesProblem.get_manufactured_solution(U, P, mu, beta, rho, t);
end

function delta = gap_size(gap, h)
if gap == "h2", delta = h^2; elseif gap == "h", delta = h; else, delta = 0; end
end

function mesh = with_tau(mesh, cfg, beta, rho_c, mu)
% "fixed": tau = 1 + rho_c B/(2 mu), B = max|beta| over the unit square;
% "class": tau = 1 (the class adds rho_c max|beta.n|/2 to S itself)
tau = 1;
if cfg == "fixed"
    [x, y] = meshgrid(linspace(0, 1, 201));
    B = max(hypot(beta{1}(x, y) + 0*x, beta{2}(x, y) + 0*x), [], 'all');
    tau = 1 + rho_c*B/(2*mu);
end
mesh = mesh.add_stabilization(tau);
end

function pb = oseen_problem(k, mu, d, cfg, y_top, varargin)
pb = IncompressibleElasticityProblem(k, "viscosity", mu, "source_term", d.f, ...
    "datum_dirichlet", d.u, "domain_measure", y_top, ...
    "convective_stabilization", cfg == "class", varargin{:});
pb.Data.convective_velocity = d.beta;
end

function [sol, mesh, global_system] = solve_strip(k, n, gap, mu, d, cfg, rho, opts, state)
% mesh of [0,1] x [0,1-delta]; points on y = 1 - delta go vertically to y = 1
delta = gap_size(gap, 1/n);
y_top = 1 - delta;
mesh = StructuredSimplicialGrid(n, n, 0, 1, 0, y_top);
mesh = with_tau(mesh, cfg, d.beta, max(rho, 1), mu);
mapping = @(X) [X(1, :); X(2, :) + (1 - X(2, :)) .* (X(2, :) > y_top - 1e-12)];
if ~isempty(state), opts = [opts, {"state", struct("fluid", state)}]; end
pb = oseen_problem(k, mu, d, cfg, y_top, "transfer_mapping", mapping, opts{:});
[sol, global_system] = solve(pb, mesh, Order = k);
end

function [sol, mesh] = run_transient(k, n, gap, mu, d_at, cfg, rho, dt)
% backward Euler from u^0 = Pi u(0); d_at(m + 1) holds the data at t = m*dt
d2 = nchoosek(k + 2, 2);
delta = gap_size(gap, 1/n);
mesh = StructuredSimplicialGrid(n, n, 0, 1, 0, 1 - delta);
u_prev = L2proj2D(mesh, d_at(1).u, k);
opts = {"density", rho, "time_step", dt, "time_scheme", "euler"};
for m = 1:numel(d_at) - 1
    [sol, mesh] = solve_strip(k, n, gap, mu, d_at(m + 1), cfg, rho, opts, struct("velocity", u_prev));
    u_prev = reshape(sol.Unknowns.field(2).dof, 2*d2, []);
end
end

function [e, r] = field_errors(sol, mesh, k, d)
% L2 errors [sigma, u, gamma, p] against the exact fields with the pressure
% shifted to zero mean over the mesh; r relative to the exact field norms
% (the velocity norm standing in for a field that vanishes)
c = domain_mean(mesh, k, d.p);
p_c = @(x, y) d.p(x, y) - c;
sigma_c = d.sigma;
sigma_c{1, 1} = @(x, y) d.sigma{1, 1}(x, y) + c;
sigma_c{2, 2} = @(x, y) d.sigma{2, 2}(x, y) + c;
e = sol.compute_l2_error({sigma_c, d.u, d.gamma});
e = [e(1:3).', pressure_error(mesh, k, sol.Unknowns.field(1).dof, p_c)];
s = [l2_norm(mesh, sigma_c(:).'), l2_norm(mesh, d.u.'), l2_norm(mesh, {d.gamma}), l2_norm(mesh, {p_c})];
s(s < 1e-12*s(2)) = s(2);
r = e ./ s;
end

function c = domain_mean(mesh, k, f)
[xy, ~, w] = quadpoints_elements_2d(mesh, k + 2, output = "cell");
c = sum(mesh.volume .* sum(w .* (f(xy{:}) + 0*xy{1}), 1)) / sum(mesh.volume);
end

function e = pressure_error(mesh, k, sigma_dof, p_ex)
% p_h = -tr(sigma_h)/2
[xy, rs, w] = quadpoints_elements_2d(mesh, k + 2, output = "array");
B = PEERSBasis(k).evaluate(rs, mesh);
S = pagemtimes(B, reshape(sigma_dof, size(B, 2), 1, 1, 1, []));
p_h = -reshape(S(:, 1, 1, 1, :) + S(:, 1, 2, 2, :), size(S, 1), [])/2;
p = reshape(p_ex(xy(1, :, :), xy(2, :, :)) + 0*xy(1, :, :), size(p_h));
e = sqrt(sum(mesh.volume .* (w'*(p_h - p).^2)));
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

function n_fail = report(name, ns, gap, E, threshold, asserted)
% table of errors and EOC; the last-pair EOC of every field is checked
% against threshold when asserted
fprintf("\n=== %s ===\n   n   delta      e(sigma)  eoc    e(u)      eoc    e(gamma)  eoc    e(p)      eoc\n", name);
eoc = [nan(1, 4); compute_eoc(E, 1./ns(:))];
for i = 1:numel(ns)
    fprintf("%4d  %.2e", ns(i), gap_size(gap, 1/ns(i)));
    fprintf("  %.3e  %5.2f", [E(i, :); eoc(i, :)]);
    fprintf("\n");
end
n_fail = 0;
if asserted && min(eoc(end, :)) < threshold
    fprintf("%s: last EOC below %.1f FAILED\n", name, threshold); n_fail = 1;
elseif ~asserted
    fprintf("%s: reported, not asserted\n", name);
end
end
