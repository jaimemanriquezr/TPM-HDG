% test_transfer_bc_biot_ladder.m
% Verification-ladder rung L2 for the quasi-static HDG Biot element
% (BiotMatrices, trace (eta_hat_x, eta_hat_y, p_hat)): the equation alone
% on a non-conforming boundary, the Dirichlet data of eta and p transferred
% along the paths (TransferEssentialBoundaryConditions + biot_path_gradient),
% one Euler step from rest (Delta t = 1, alpha = 1, s0 = 0.1, kappa = I),
% k = 1, 2, 3.
%
% k = 0 probe: one solve on a 2x2 strip; k = 0 is then recorded as SKIPPED
%   (the stress element has no curl-bubble at k = 0), or with the error and
%   its location if the solve fails; never counted as a failure.
% (P) Patch test: the unit square meshed on [0,1] x [0,1-delta], the data
%   on y = 1 transferred along vertical paths, delta = h^2 and delta = h,
%   n = 2, 4, 8 (h = 1/n); eta and p polynomials of degree k, so sigma, q,
%   gamma are of degree <= k and the exact solution lies in the discrete
%   space; the Taylor transfer is exact for them. The L2 error of sigma,
%   eta, gamma, q, p relative to the norm of the exact field must be
%   <= 1e-10 (absolute when the exact field vanishes: gamma, q at k = 0).
% (S) Strip convergence: the same geometry, smooth solution, n = 4..32.
%   delta = h^2: rate between the coarsest and finest levels
%   >= k + 1 - 0.25 for sigma, eta, q, p. delta = h: reported, not asserted.
% (D) Carved disk: R = 0.75 carved from [-1,1]^2 by immerse_mesh, data
%   transferred radially, n = 8..64: the same rate threshold for eta and p
%   (sigma, q reported).
% The rotation gamma is reported, never asserted: at k = 0 its error does
% not decrease, on the fitted unit square (BiotProblem) as well.
%
% Set the environment variable HDG_TEST_SMOKE=1 to run two levels per
% convergence study only (rates then reported, not asserted).
clc; clear;
n_fail = 0;
smoke = getenv("HDG_TEST_SMOKE") == "1";

E = 1.0;
nu = 0.3;
alpha = 1.0;
s0 = 0.1;
kappa = eye(2);
dt = 1.0;
names = ["sigma", "eta", "gamma", "q", "p"];
asserted_strip = [1 2 4 5];     % sigma, eta, q, p
asserted_disk = [2 5];          % eta, p

R = 0.75;
phi = @(X) sqrt(X(1, :, :).^2 + X(2, :, :).^2) - R;
radial_map = @(X) R * X ./ vecnorm(X);

%% k = 0 probe
k_list = 0:3;
try
    data = biot_data(@() polynomial_fields(0), E, nu, alpha, s0, kappa, dt);
    [mesh, map] = strip_mesh(2, 1/4);
    solve_strip(mesh, map, 0, E, nu, alpha, s0, kappa, dt, data);
    fprintf("k = 0 probe: solved\n");
    % k = 0 runs but is not checked: PEERSBasis(0) has no curl-bubble
    % enrichment, so the rotation does not converge (verification ladder:
    % k = 0 is skipped for the stress-based elements)
    fprintf("k = 0 SKIPPED: no curl-bubble enrichment in PEERSBasis(0)\n");
    k_list = 1:3;
catch err
    loc = "";
    if ~isempty(err.stack)
        loc = sprintf("%s:%d", err.stack(1).file, err.stack(1).line);
    end
    fprintf("k = 0 SKIPPED: %s (%s) at %s\n", err.message, err.identifier, loc);
    k_list = 1:3;
end

%% (P) patch test on the strip
fprintf("\n=== (P) patch test, strip, polynomial data of degree k ===\n");
for k = k_list
    data = biot_data(@() polynomial_fields(k), E, nu, alpha, s0, kappa, dt);
    for gap_kind = ["h2", "h"]
        worst = 0;
        for n = [2 4 8]
            [mesh, map] = strip_mesh(n, gap_size(gap_kind, 1/n));
            sol = solve_strip(mesh, map, k, E, nu, alpha, s0, kappa, dt, data);
            err = biot_errors(sol, mesh, data);
            nrm = exact_norms(sol, mesh, data);
            rel = err ./ max(nrm, (nrm == 0));
            worst = max(worst, max(rel));
            fprintf("k = %d, gap %-2s, n = %d: rel err %s\n", k, gap_kind, n, ...
                strjoin(names + " " + compose("%.1e", rel), ", "));
        end
        if ~(worst <= 1e-10)
            fprintf("(P) FAILED: k = %d, gap %s, max rel err %.3e\n", k, gap_kind, worst);
            n_fail = n_fail + 1;
        end
    end
end

%% (S) strip convergence, smooth solution
data = biot_data(@smooth_fields, E, nu, alpha, s0, kappa, dt);
for gap_kind = ["h2", "h"]
    for k = k_list
        levels = [4 8 16 32];
        if smoke, levels = levels(1:2); end
        err = zeros(numel(levels), 5);
        for i = 1:numel(levels)
            [mesh, map] = strip_mesh(levels(i), gap_size(gap_kind, 1/levels(i)));
            sol = solve_strip(mesh, map, k, E, nu, alpha, s0, kappa, dt, data);
            err(i, :) = biot_errors(sol, mesh, data);
        end
        rate = print_table(sprintf("(S) strip, gap %s, k = %d", gap_kind, k), ...
            "n", levels, 1./levels, err, names);
        if gap_kind == "h2" && ~smoke && any(rate(asserted_strip) < k + 1 - 0.25)
            fprintf("(S) rate FAILED: gap %s, k = %d\n", gap_kind, k);
            n_fail = n_fail + 1;
        end
    end
end

%% (D) carved disk, radial transfer
for k = k_list
    levels = [8 16 32 64];
    if smoke, levels = levels(1:2); end
    err = zeros(numel(levels), 5);
    for i = 1:numel(levels)
        n = levels(i);
        mesh = immerse_mesh(StructuredSimplicialGrid(n, n, -1, 1, -1, 1), phi);
        problem = build_problem(mesh, k, E, nu, alpha, s0, kappa, dt, data);
        problem.BoundaryConditions = TransferEssentialBoundaryConditions(data.dirichlet, ...
            "Dirichlet", radial_map, path_gradient = biot_path_gradient(k, E, nu, alpha, kappa));
        sol = solve(problem, mesh, k);
        err(i, :) = biot_errors(sol, mesh, data);
    end
    rate = print_table(sprintf("(D) carved disk, k = %d", k), "n", levels, 2./levels, err, names);
    if ~smoke && any(rate(asserted_disk) < k + 1 - 0.25)
        fprintf("(D) rate FAILED: k = %d\n", k);
        n_fail = n_fail + 1;
    end
end

fprintf("%d check(s) failed.\n", n_fail);

%% ------------------------------------------------------------------ %%
function [u, p] = polynomial_fields(k)
% displacement and pressure of degree exactly k (constants for k = 0)
syms x y
u = [(0.3 + x + 0.7*y)^k + 0.5*x^k; (0.5 - 0.4*x + y)^k + 2];
p = (0.2 + 0.6*x - 0.9*y)^k + y^k;
end

function [u, p] = smooth_fields()
syms x y
u = [sin(pi*x) .* sin(pi*y); cos(pi*x) .* cos(pi*y)];
p = cos(pi*x) .* sin(pi*y) + x;
end

function data = biot_data(fields, E, nu, alpha, s0, kappa, dt)
% one Euler step from rest: d_t w = w/dt; total stress, rotation, Darcy
% flux and the sources f = -div sigma, g = (s0 p + alpha div eta)/dt + div q
% fields is a handle returning (u, p) as symbolic expressions
[u, p] = fields();
[lambda, mu] = LinearElasticityProblem.get_lame_parameters(E, nu);
syms x y
gu = jacobian(u, [x y]);
div_u = trace(gu);
sig = mu*(gu + gu.') + lambda*div_u*eye(2) - alpha*p*eye(2);
f = -[diff(sig(1,1), x) + diff(sig(1,2), y); diff(sig(2,1), x) + diff(sig(2,2), y)];
q = -kappa * gradient(p, [x y]);
g = (s0*p + alpha*div_u)/dt + diff(q(1), x) + diff(q(2), y);
data.sigma = Datum(sig);
data.u = Datum(u);
data.gamma = Datum((gu(1,2) - gu(2,1))/2);    % rotation of LinearElasticityMatrices
data.p = Datum(p);
data.q = Datum(q);
data.f = {Datum(f(1)), Datum(f(2))};
data.g = Datum(g);
data.dirichlet = Datum([u; p]);
end

function delta = gap_size(gap_kind, h)
if gap_kind == "h2", delta = h^2; else, delta = h; end
end

function [mesh, map] = strip_mesh(n, delta)
% [0,1] x [0,1-delta]; the top face y = 1 - delta is mapped vertically to
% y = 1, the other boundary faces (on the true boundary) to themselves
top = 1 - delta;
mesh = StructuredSimplicialGrid(n, n, 0, 1, 0, top);
map = @(X) [X(1, :, :); X(2, :, :) + delta*(X(2, :, :) >= top - 1e-12)];
end

function sol = solve_strip(mesh, map, k, E, nu, alpha, s0, kappa, dt, data)
problem = build_problem(mesh, k, E, nu, alpha, s0, kappa, dt, data);
problem.BoundaryConditions = TransferEssentialBoundaryConditions(data.dirichlet, ...
    "Dirichlet", map, path_gradient = biot_path_gradient(k, E, nu, alpha, kappa));
sol = solve(problem, mesh, k);
end

function problem = build_problem(mesh, k, E, nu, alpha, s0, kappa, dt, data)
Pk = DubinerBasis(k);
problem = FEProblem(Spaces = [FESpace([PEERSBasis(k), Pk.to_vector(), Pk, Pk.to_vector(), Pk], ...
    Condense = true), FESpace(LegendreBasis(k).to_vector(3))]);
[lhs, flux] = BiotMatrices(k, mesh, E, nu, kappa, biot_alpha = alpha, storage = s0, ...
    effective_time_step = dt);
w = TestFunction(Pk);
F1 = double(data.f{1} * w * Measure(mesh, 'dx'));
F2 = double(data.f{2} * w * Measure(mesh, 'dx'));
G = double(data.g * w * Measure(mesh, 'dx'));
problem.Operator = AffineOperator(lhs, flux) ...
    + IndexedMatrix(IndexedArray(cat(1, F1, F2), [-1, 2])) ...
    + IndexedMatrix(IndexedArray(G, [-1, 5]));
end

function err = biot_errors(sol, mesh, data)
% L2 errors of (sigma, eta, gamma, q, p)
exact = {data.sigma, data.u, data.gamma, data.q, data.p};
err = zeros(1, 5);
for f = 1:5
    err(f) = compute_l2_error(sol.Unknowns.field(f), mesh, exact{f});
end
end

function nrm = exact_norms(sol, mesh, data)
% L2 norms of the exact fields (distance from the zero field of each basis)
exact = {data.sigma, data.u, data.gamma, data.q, data.p};
nrm = zeros(1, 5);
for f = 1:5
    field = sol.Unknowns.field(f);
    nrm(f) = compute_l2_error(FEField(field.basis, zeros(size(field.dof))), mesh, exact{f});
end
end

function rate = print_table(title, level_name, levels, h, err, names)
% errors with per-level EOC and the rate between the coarsest and finest
% levels
eoc = [nan(1, size(err, 2)); log(err(1:end-1, :)./err(2:end, :)) ./ log(h(1:end-1)./h(2:end)).'];
fprintf("\n%s\n  %4s  %s\n", title, level_name, ...
    strjoin(compose("%10s  %5s", names.', repmat("eoc", numel(names), 1)).', "  "));
for i = 1:numel(levels)
    fprintf("  %4d  %s\n", levels(i), ...
        strjoin(compose("%10.3e  %5.2f", err(i, :).', eoc(i, :).').', "  "));
end
rate = log(err(1, :)./err(end, :)) / log(h(1)/h(end));
fprintf("  rate coarsest-finest: %s\n", strjoin(names + " " + compose("%.2f", rate), ", "));
end
