% test_incompressible_elasticity_transient_ladder.m
% Verification ladder, rung L1, for the TRANSIENT stress-form Stokes fluid
% (PEERS stress sigma, velocity u, rotation gamma; density rho > 0, no
% convection) alone on the unit square, k = 1, 2, 3, on two assembly paths
% under the "pinned" closure (traction face x = 1 with datum sigma n):
%   matrix: FEProblem + IncompressibleElasticityMatrices + AffineOperator
%           with ComponentwiseBoundaryConditions, the fluid of
%           StokesBiotCoupling.interface_problem;
%   class:  solve(IncompressibleElasticityProblem, mesh, k).
% Every step is one backward-Euler-type solve
%   rho (u^{n+1} - w*)/cdt - div sigma^{n+1} = f(t^{n+1}),
% passed as density = rho, theta = 1, time_step = cdt, state.velocity = w*
% (class: state.fluid.velocity = w*), with
%   euler:        cdt = dt,      w* = u^n;
%   BDF2 slot:    cdt = 2 dt/3,  w* = (4 u^n - u^{n-1})/3,
% so that (u^{n+1} - w*)/cdt = (3 u^{n+1} - 4 u^n + u^{n-1})/(2 dt). The
% BDF2 slot starts either from the exact two-level start u^0 = Pi u(0),
% u^1 = Pi u(dt) ("bdf2-exact") or from one Euler step ("bdf2-euler").
% Pi is the L2 projection; f contains rho u_t(t^{n+1}) exactly.
% p_h = -tr(sigma_h)/2 is reported as its own field.
%
%  1. Patch test: velocity and pressure polynomials of degree <= k in
%     space (stream function of degree k + 1; sigma, gamma, f of degree
%     <= k). Euler: u linear in t, 3 steps; the backward difference is
%     then exact. BDF2 slot: u quadratic in t from the exact two-level
%     start, 3 steps; the three-point backward difference is exact for
%     quadratics. The exact solution then solves every discrete step: max
%     over steps of the relative L2 error of sigma, u, gamma, p <= 1e-10
%     on three meshes, both paths.
%  2. Space ladder, smooth solution, n = 4, 8, 16, 32, dt = 0.1: euler
%     with data linear in t (3 steps) and the BDF2 slot with data
%     quadratic in t from the exact start (3 steps after the start), so
%     the errors are spatial only; the EOC of the last pair at the final
%     time must be >= k + 1 - 0.3 for sigma, u, gamma and p.
%  3. Time order (k = 2, 4 x 4 mesh, matrix path): u and p of degree <= k
%     in space and smooth in t (sin(3t + 1), cos(2t + 1): u_tt(0) is not
%     zero, so the Euler start of bdf2-euler is not artificially
%     accurate), T = 1, dt = 1/4 ... 1/64. The
%     errors are then temporal (the spatial part is O(h^{k+1}) times the
%     temporal one). Errors at t = T and the max over the steps, all
%     fields. Asserted, EOC of the last pair at t = T, all fields: euler
%     >= 0.9, bdf2-exact >= 1.8. bdf2-euler is reported only.
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
    patch_meshes = {[2 2], [3 5], [6 4]};
    for k = ks
        for scheme = ["euler", "bdf2-exact"]
            tU = time_U.(extractBefore(scheme + "-", "-"));
            n_levels = 4 + (scheme ~= "euler");
            D = data_at(tU + "*(" + patch_U(k + 1) + ")", time_P + "*(" + patch_P(k + 1) + ")", ...
                mu, rho, (0:n_levels - 1)*dt);
            for path = paths
                worst = 0;
                for m = 1:numel(patch_meshes)
                    nm = patch_meshes{m};
                    mesh = unit_square(nm, tau);
                    [E, S] = march(path, k, mesh, D, mu, rho, dt, scheme);
                    worst = max(worst, max(E ./ S, [], "all", "omitnan"));
                end
                fprintf("patch k = %d %-10s %-6s: max rel. error over fields, steps, meshes %.2e\n", ...
                    k, scheme, path, worst);
                if ~(worst <= 1e-10)
                    fprintf("patch k = %d %s %s FAILED\n", k, scheme, path);
                    n_fail = n_fail + 1;
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
            for path = paths
                E = zeros(numel(ns), 4);
                for i = 1:numel(ns)
                    mesh = unit_square([ns(i) ns(i)], tau);
                    Ei = march(path, k, mesh, D, mu, rho, dt, scheme);
                    E(i, :) = Ei(end, :);
                end
                name = sprintf("space %s %s k = %d (t = %.1f)", scheme, path, k, (n_levels - 1)*dt);
                n_fail = n_fail + report(name, "N", ns, 1./ns, E, k + 1 - 0.3, true(1, 4));
            end
        end
    end
end

%% 3. time order
if any(parts == "time") && any(ks == 2)
    k = 2;
    mesh = unit_square([4 4], tau);
    dts = 1./[4 8 16 32 64];
    D_fine = data_at("(2 + sin(3*t + 1))*(" + patch_U(k + 1) + ")", "cos(2*t + 1)*(" + patch_P(k + 1) + ")", ...
        mu, rho, (0:64)/64);
    expected = struct("euler", 0.9, "bdf2_exact", 1.8, "bdf2_euler", NaN);
    for scheme = schemes
        [E_T, E_max, E_1] = deal(zeros(numel(dts), 4));
        for i = 1:numel(dts)
            D = D_fine(1:round(64*dts(i)):end);
            E = march("matrix", k, mesh, D, mu, rho, dts(i), scheme);
            E_T(i, :) = E(end, :);
            E_max(i, :) = max(E, [], 1, "omitnan");
            E_1(i, :) = E(find(~isnan(E(:, 1)), 1), :);
        end
        target = expected.(strrep(scheme, "-", "_"));
        asserted = repmat(~isnan(target), 1, 4);
        n_fail = n_fail + report(sprintf("time %s k = %d, errors at t = T = 1", scheme, k), ...
            "1/dt", 1./dts, dts, E_T, target, asserted);
        report(sprintf("time %s k = %d, max over the steps (reported)", scheme, k), ...
            "1/dt", 1./dts, dts, E_max, NaN, false(1, 4));
        report(sprintf("time %s k = %d, errors after the first computed step (reported)", scheme, k), ...
            "1/dt", 1./dts, dts, E_1, NaN, false(1, 4));
    end
end

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_incompressible_elasticity_transient_ladder:failed", "%d checks failed", n_fail);
end

%% ------------------------------------------------------------------------
function D = data_at(U, P, mu, rho, ts)
D = arrayfun(@(t) IncompressibleElasticityProblem.get_manufactured_solution(U, P, mu, rho, t), ts);
end

function mesh = unit_square(nm, tau)
mesh = StructuredSimplicialGrid(nm(1), nm(2)).add_stabilization(tau);
mesh = apply_tag(mesh, "Neumann", @(X) isapprox(X(1, :), 1));
end

function [E, S] = march(path, k, mesh, D, mu, rho, dt, scheme)
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
    sol = solve_step(path, k, mesh, D(n + 1), mu, rho, cdt, w);
    E(n, :) = field_errors(sol, mesh, k, D(n + 1));
    if nargout > 1, S(n, :) = field_scales(mesh, D(n + 1)); end
    u_old = u_now;
    u_now = reshape(sol.Unknowns.field(2).dof, 2*d2, []);
end
end

function sol = solve_step(path, k, mesh, d, mu, rho, cdt, w)
% one step rho (u - w)/cdt - div sigma = f with Dirichlet datum u on the
% "Dirichlet" faces and traction sigma n on x = 1, closure "pinned"
gN = {d.sigma{1,1}; d.sigma{2,1}};
if path == "matrix"
    Pk = DubinerBasis(k);
    pb = FEProblem(Order = k, Spaces = [ ...
        FESpace([PEERSBasis(k), Pk.to_vector(), Pk, DubinerBasis(0)], Condense = true), ...
        FESpace([LegendreBasis(k).to_vector(2), DubinerBasis(0)])]);
    pb.Mesh = mesh;
    bc = ComponentwiseBoundaryConditions([ ...
        boundary_condition("Dirichlet", x_data = d.u{1}, y_data = d.u{2});
        boundary_condition("Neumann", x = "neumann", y = "neumann", ...
            x_data = gN{1}, y_data = gN{2})]);
    [lhs, flux, rhs] = IncompressibleElasticityMatrices(k, mesh, mu, d.f, ...
        "pressure_mean", "pinned", "natural_faces", bc.natural_faces(mesh), ...
        "density", rho, "theta", 1, "time_step", cdt, "state", struct("velocity", w));
    pb.Operator = AffineOperator(lhs + rhs) + flux;
    pb.BoundaryConditions = bc;
    sol = solve(pb, mesh, k);
else
    pb = IncompressibleElasticityProblem(k, viscosity = mu, pressure_mean = "pinned", ...
        source_term = d.f, datum_dirichlet = d.u, datum_neumann = gN, ...
        density = rho, theta = 1, time_step = cdt, state = struct("fluid", struct("velocity", w)));
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
