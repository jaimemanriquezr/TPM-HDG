% test_transfer_bc_incompressible_elasticity_transient.m
% Verification ladder, rung L2, for the TRANSIENT stress-form Stokes fluid
% (density rho > 0, no convection) alone through TPM, k = 1, 2, 3: the
% domain is the unit square, the mesh covers the strip [0,1] x [0,1-delta]
% (delta = h^2 and delta = h, h = 1/N), and the velocity datum on y = 1 is
% transferred along vertical paths of length delta to the top face (the
% paths of the other "Dirichlet" faces have length zero); the face x = 1
% carries the traction sigma n (closure "pinned"). Two assembly paths:
%   matrix: FEProblem + IncompressibleElasticityMatrices + AffineOperator
%           with [ComponentwiseBoundaryConditions (traction),
%           TransferEssentialBoundaryConditions.incompressible_stress
%           (Dirichlet)], the fluid of StokesBiotCoupling.interface_problem
%           with the transfer of test_transfer_bc_incompressible_elasticity;
%   class:  IncompressibleElasticityProblem with transfer_mapping
%           (apply_transfer_conditions).
% Time stepping and the three schemes (euler, bdf2-exact, bdf2-euler) as
% in test_incompressible_elasticity_transient_ladder: one solve per step
% with density = rho, theta = 1, time_step = cdt, state.velocity = w*,
%   euler: cdt = dt, w* = u^n;  BDF2 slot: cdt = 2 dt/3,
%   w* = (4 u^n - u^{n-1})/3.
% The history is volume data of the strip mesh at t_n, the transferred
% datum is boundary data on y = 1 at t_{n+1}. p_h = -tr(sigma_h)/2.
%
%  1. Patch test: u, p polynomials of degree <= k in space (stream
%     function of degree k + 1), u linear in t (euler, 3 steps) or
%     quadratic in t (bdf2-exact, 3 steps after the exact start); the
%     backward differences and the Taylor transfer are then exact: max
%     over steps of the relative L2 error of sigma, u, gamma, p <= 1e-10,
%     N = 4, 8, 16, both gaps, both paths. Not asserted for k = 3 with
%     delta = h: there the homogeneous step (zero data) has growing modes,
%     |g| about 2.2-2.5 per Euler step of dt = 0.1 and 1.2-3.6 at dt =
%     0.01 on N = 4, 8, 16 (power iteration, both paths alike), so the
%     round-off of the exact solution grows by about 2x per Euler step and
%     20-35x per BDF2 step (patch error 6e-10 / 7e-8 after 3 steps, 1e-8 /
%     1e-1 after 8). The steady solve on the same meshes is exact to 1e-11:
%     the transferred rows make the discrete operator non-dissipative for
%     paths as long as an element. With delta = h^2 every mode decays.
%  2. Space ladder, smooth solution, N = 4, 8, 16, 32, dt = 0.1, euler
%     (linear in t) and bdf2-exact (quadratic in t): errors at the final
%     time; delta = h^2: last-pair EOC >= k + 1 - 0.3 for sigma and u (the
%     thresholds of test_transfer_bc_incompressible_elasticity; gamma and
%     p printed); delta = h: reported only (with k = 3 the BDF2 errors are
%     erratic in N, a consequence of the growing modes above).
%  3. Time order (k = 2, N = 4, delta = h^2, matrix path): u, p of degree
%     <= k in space, smooth in t (sin(3t + 1), cos(2t + 1); u_tt(0) is not
%     zero), T = 1, dt = 1/4 ... 1/64; errors at t = T
%     (asserted, last-pair EOC, all fields: euler >= 0.9, bdf2-exact >=
%     1.8; bdf2-euler reported), the max over the steps and the first
%     computed step (reported).
%
% Pieces: HDG_LADDER_K (e.g. "2" or "1,3") restricts k; HDG_LADDER_PARTS
% (any of "patch,space,time") restricts the parts. Default: everything.
clc; clear;
n_fail = 0;

mu = 0.1; tau = 1; rho = 1;
ks = env_numbers("HDG_LADDER_K", 1:3);
parts = env_strings("HDG_LADDER_PARTS", ["patch", "space", "time"]);
paths = ["matrix", "class"];
schemes = ["euler", "bdf2-exact", "bdf2-euler"];
patch_U = ["2*x - 3*y", ...
           "x^2 + 3*x*y - y^2 + 2*x", ...
           "x^3 + 2*x^2*y - x*y^2 + y^3 + x*y", ...
           "x^4 - 2*x^3*y + x^2*y^2 + y^4 + x^2*y"];
patch_P = ["0", "x - 2*y", "x^2 - x*y + 3*y", "x^3 - x*y^2 + y^2 + x"];
% time factors: linear (euler exact) and quadratic (BDF2 exact)
time_U = struct("euler", "(1 + t)", "bdf2", "(1 + t - 2*t^2)");
time_P = "(1 + t^2)";

%% 1. patch test
if any(parts == "patch")
    dt = 0.1;
    for k = ks
        for scheme = ["euler", "bdf2-exact"]
            tU = time_U.(extractBefore(scheme + "-", "-"));
            n_levels = 4 + (scheme ~= "euler");
            D = data_at(tU + "*(" + patch_U(k + 1) + ")", time_P + "*(" + patch_P(k + 1) + ")", ...
                mu, rho, (0:n_levels - 1)*dt);
            for gap = ["h2", "h"]
                for path = paths
                    worst = 0;
                    for N = [4 8 16]
                        [mesh, mapping] = strip_mesh(N, gap, tau);
                        [E, S] = march(path, k, mesh, mapping, D, mu, rho, dt, scheme);
                        worst = max(worst, max(E ./ S, [], "all", "omitnan"));
                    end
                    fprintf("patch k = %d %-10s gap %-2s %-6s: max rel. error over fields, steps, meshes %.2e\n", ...
                        k, scheme, gap, path, worst);
                    if k == 3 && gap == "h"
                        % reported only: the step operator has growing modes (header)
                        fprintf("patch k = %d %s gap h %s: NOT ASSERTED (growing modes, see header)\n", ...
                            k, scheme, path);
                    elseif ~(worst <= 1e-10)
                        fprintf("patch k = %d %s gap %s %s FAILED\n", k, scheme, gap, path);
                        n_fail = n_fail + 1;
                    end
                end
            end
        end
    end
end

%% 2. space ladder
if any(parts == "space")
    dt = 0.1;
    ns = [4 8 16 32];
    smooth_U = "sin(pi*x + 0.3)*sin(pi*y + 0.7)";
    smooth_P = "cos(pi*x)*sin(pi*y) + x*y";
    for scheme = ["euler", "bdf2-exact"]
        tU = time_U.(extractBefore(scheme + "-", "-"));
        n_levels = 4 + (scheme ~= "euler");
        D = data_at(tU + "*" + smooth_U, time_P + "*(" + smooth_P + ")", mu, rho, (0:n_levels - 1)*dt);
        for k = ks
            for gap = ["h2", "h"]
                for path = paths
                    E = zeros(numel(ns), 4);
                    for i = 1:numel(ns)
                        [mesh, mapping] = strip_mesh(ns(i), gap, tau);
                        Ei = march(path, k, mesh, mapping, D, mu, rho, dt, scheme);
                        E(i, :) = Ei(end, :);
                    end
                    name = sprintf("space %s gap %s %s k = %d (t = %.1f)%s", scheme, gap, path, k, ...
                        (n_levels - 1)*dt, ternary(gap == "h", " [reported only]", ""));
                    % gap h^2: sigma and u asserted (gamma, p printed); gap h: reported
                    n_fail = n_fail + report(name, "N", ns, 1./ns, E, k + 1 - 0.3, ...
                        (gap == "h2") & [true, true, false, false]);
                end
            end
        end
    end
end

%% 3. time order
if any(parts == "time") && any(ks == 2)
    k = 2;
    [mesh, mapping] = strip_mesh(4, "h2", tau);
    dts = 1./[4 8 16 32 64];
    D_fine = data_at("(2 + sin(3*t + 1))*(" + patch_U(k + 1) + ")", "cos(2*t + 1)*(" + patch_P(k + 1) + ")", ...
        mu, rho, (0:64)/64);
    expected = struct("euler", 0.9, "bdf2_exact", 1.8, "bdf2_euler", NaN);
    for scheme = schemes
        [E_T, E_max, E_1] = deal(zeros(numel(dts), 4));
        for i = 1:numel(dts)
            D = D_fine(1:round(64*dts(i)):end);
            E = march("matrix", k, mesh, mapping, D, mu, rho, dts(i), scheme);
            E_T(i, :) = E(end, :);
            E_max(i, :) = max(E, [], 1, "omitnan");
            E_1(i, :) = E(find(~isnan(E(:, 1)), 1), :);
        end
        target = expected.(strrep(scheme, "-", "_"));
        asserted = repmat(~isnan(target), 1, 4);
        n_fail = n_fail + report(sprintf("time %s k = %d gap h2 N = 4, errors at t = T = 1", scheme, k), ...
            "1/dt", 1./dts, dts, E_T, target, asserted);
        report(sprintf("time %s k = %d gap h2 N = 4, max over the steps (reported)", scheme, k), ...
            "1/dt", 1./dts, dts, E_max, NaN, false(1, 4));
        report(sprintf("time %s k = %d gap h2 N = 4, errors after the first computed step (reported)", scheme, k), ...
            "1/dt", 1./dts, dts, E_1, NaN, false(1, 4));
    end
end

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_transfer_bc_incompressible_elasticity_transient:failed", "%d checks failed", n_fail);
end

%% ------------------------------------------------------------------------
function D = data_at(U, P, mu, rho, ts)
D = arrayfun(@(t) IncompressibleElasticityProblem.get_manufactured_solution(U, P, mu, rho, t), ts);
end

function [mesh, mapping] = strip_mesh(N, gap, tau)
% mesh of [0,1] x [0,1-delta], "Neumann" on x = 1; points on y = 1 - delta
% go vertically to y = 1, the others stay
if gap == "h2", delta = 1/N^2; else, delta = 1/N; end
y_top = 1 - delta;
mesh = StructuredSimplicialGrid(N, N, 0, 1, 0, y_top).add_stabilization(tau);
mesh = apply_tag(mesh, "Neumann", @(X) isapprox(X(1, :), 1));
mapping = @(X) [X(1, :); X(2, :) + (1 - X(2, :)) .* (X(2, :) > y_top - 1e-12)];
end

function [E, S] = march(path, k, mesh, mapping, D, mu, rho, dt, scheme)
% D(j + 1) holds the data at t = j*dt. Rows of E: L2 errors [sigma, u,
% gamma, p] after step n = 1 .. numel(D) - 1 (NaN for the exact start
% level u^1 of bdf2-exact); S: the norms of the exact fields
d2 = nchoosek(k + 2, 2);
n_steps = numel(D) - 1;
[E, S] = deal(NaN(n_steps, 4));
u_old = [];
u_now = L2proj2D(mesh, D(1).u, k);
n_first = 1;
if scheme == "bdf2-exact"
    u_old = u_now;
    u_now = L2proj2D(mesh, D(2).u, k);
    n_first = 2;
end
for n = n_first:n_steps
    if scheme == "euler" || isempty(u_old)
        cdt = dt;  w = u_now;
    else
        cdt = 2*dt/3;  w = (4*reshape(u_now, 2*d2, []) - reshape(u_old, 2*d2, []))/3;
    end
    sol = solve_step(path, k, mesh, mapping, D(n + 1), mu, rho, cdt, w);
    E(n, :) = field_errors(sol, mesh, k, D(n + 1));
    if nargout > 1, S(n, :) = field_scales(mesh, D(n + 1)); end
    u_old = u_now;
    u_now = reshape(sol.Unknowns.field(2).dof, 2*d2, []);
end
end

function sol = solve_step(path, k, mesh, mapping, d, mu, rho, cdt, w)
% one step rho (u - w)/cdt - div sigma = f with the Dirichlet datum u
% transferred to the "Dirichlet" faces and traction sigma n on x = 1,
% closure "pinned"
gN = {d.sigma{1,1}; d.sigma{2,1}};
if path == "matrix"
    Pk = DubinerBasis(k);
    pb = FEProblem(Order = k, Spaces = [ ...
        FESpace([PEERSBasis(k), Pk.to_vector(), Pk, DubinerBasis(0)], Condense = true), ...
        FESpace([LegendreBasis(k).to_vector(2), DubinerBasis(0)])]);
    pb.Mesh = mesh;
    traction = ComponentwiseBoundaryConditions(boundary_condition("Neumann", ...
        x = "neumann", y = "neumann", x_data = gN{1}, y_data = gN{2}));
    bc = [traction, TransferEssentialBoundaryConditions.incompressible_stress( ...
        Datum(d.u), "Dirichlet", mapping, mu)];
    [lhs, flux, rhs] = IncompressibleElasticityMatrices(k, mesh, mu, d.f, ...
        "pressure_mean", "pinned", "natural_faces", traction.natural_faces(mesh), ...
        "density", rho, "theta", 1, "time_step", cdt, "state", struct("velocity", w));
    pb.Operator = AffineOperator(lhs + rhs) + flux;
    pb.BoundaryConditions = bc;
    sol = solve(pb, mesh, k);
else
    pb = IncompressibleElasticityProblem(k, viscosity = mu, pressure_mean = "pinned", ...
        source_term = d.f, datum_dirichlet = d.u, datum_neumann = gN, ...
        density = rho, theta = 1, time_step = cdt, state = struct("fluid", struct("velocity", w)), ...
        transfer_mapping = mapping);
    sol = solve(pb, mesh, k);
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

function n_fail = report(name, label, xs, hs, E, threshold, asserted)
% EOC table of the errors E (rows over hs); the last-pair EOC of the
% asserted fields must be >= threshold
fprintf("%s\n  %4s   |sigma|     EOC    |u|         EOC    |gamma|     EOC    |p|         EOC\n", name, label);
eoc = compute_eoc(E, hs(:));
r = NaN(1, 4);
for i = 1:numel(xs)
    r = NaN(1, 4);
    if i > 1, r = eoc(i - 1, :); end
    fprintf("  %4d   %.3e  %5.2f   %.3e  %5.2f   %.3e  %5.2f   %.3e  %5.2f\n", xs(i), ...
        E(i, 1), r(1), E(i, 2), r(2), E(i, 3), r(3), E(i, 4), r(4));
end
n_fail = double(any(asserted) && any(r(asserted) < threshold));
if n_fail
    fprintf("%s: last EOC below %.1f FAILED\n", name, threshold);
end
end

function v = ternary(c, a, b)
if c, v = a; else, v = b; end
end

function v = env_numbers(name, default)
v = default;
s = getenv(name);
if s ~= "", v = str2double(split(string(s), ",")).'; end
end

function v = env_strings(name, default)
v = default;
s = getenv(name);
if s ~= "", v = strtrim(split(string(s), ",")).'; end
end
