% test_initial_pressure_fluid.m
% Initial state of the stress-velocity fluid (IncompressibleElasticity /
% NavierStokesProblem.initial_state): the acceleration problem
%   rho a + grad p = f(0) + div sigma_v [- rho (beta.grad) u],  div a = 0,
% solved for (u_t(0), p(0), sigma(0)) with the velocity-hybridized HDG of
% AccelerationMatrices.
%  1. Rates of u_t0, p0, sigma0 against the P_k projections of the
%     manufactured u_t(0), p(0), sigma(0) (relative L2), fitted unit square
%     (Stokes and Navier-Stokes; constrained = all Dirichlet, pinned = a
%     traction face x = 1 with the nonzero datum sigma(0) n), and the TPM
%     transfer on an immersed disk (radial mapping). Asserted at the
%     measured order minus 0.3 on the last refinement; measured orders:
%       fitted Stokes, both closures       u_t0, p0, sigma0: k+1
%       fitted Navier-Stokes, constrained  u_t0, p0, sigma0: k+1
%       fitted Navier-Stokes, pinned       u_t0: k+1/2 (no convective flux
%                                          on the traction face), p0, sigma0: k+1
%       transfer (disk), Stokes            k = 1: 2;  k = 2: 2 (not k+1)
%  2. Polynomial reproduction (u, p of degree <= k): u_t0, p0, sigma0 to
%     round-off, fitted (both closures, Stokes k <= 3 and Navier-Stokes
%     k <= 2) and through the transfer (k <= 2).
%  3. theta = 1/2: one and four Crank-Nicolson steps from initial_state vs
%     from the exact Pi u_t(0); the difference shrinks with h.
%  4. Laziness: theta = 1 builds no acceleration problem unless
%     with_pressure = true.
%  5. Round trip: initial_state -> step -> get_kinematic_state equals the
%     same step from a hand-built state, bitwise; the pressure of the advanced state is -tr(sigma)/2.
clc; clear;
n_fail = 0;
mu = 0.1; rho = 1.3; dt = 1e-2;
U = "(x^2*y^2/2 + sin(pi*x)*sin(pi*y)/pi)*cos(t + x*y)";
P = "sin(pi*x)*cos(pi*y)*(1 + t) + x^2*exp(t)";
dS = IncompressibleElasticityProblem.get_manufactured_solution(U, P, mu, rho, 0);
dN = NavierStokesProblem.get_manufactured_solution(U, P, mu, [], rho, 0);
R = 0.75;
phi = @(X) sqrt(X(1, :, :).^2 + X(2, :, :).^2) - R;
radial = @(X) R * X ./ vecnorm(X);

%% 1. rates
for arm = ["stokes", "ns"]
    for closure = ["constrained", "pinned"]
        for k = 1:3
            ns = [4 8 16 32];
            if k == 3, ns = [4 8 16]; end
            E = zeros(numel(ns), 3);
            for i = 1:numel(ns)
                mesh = StructuredSimplicialGrid(ns(i), ns(i)).add_stabilization(1);
                d = ternary(arm == "stokes", dS, dN);
                [pde, mesh] = fluid_problem(arm, closure, k, mesh, d, mu, rho, dt, "cranknicolson");
                st = initial_state(pde, mesh, k);
                E(i, :) = fluid_errors(mesh, k, st.fluid, d, false);
            end
            observed = [k + 1, k + 1, k + 1];
            if arm == "ns" && closure == "pinned", observed(1) = k + 1/2; end
            n_fail = n_fail + report(sprintf("%s %s k = %d", arm, closure, k), ns, E, observed);
        end
    end
end
for k = 1:2
    ns = [8 16 32];
    E = zeros(numel(ns), 3);
    for i = 1:numel(ns)
        mesh = immerse_mesh(StructuredSimplicialGrid(ns(i), ns(i), -1, 1, -1, 1), phi).add_stabilization(1);
        [pde, mesh] = fluid_problem("stokes", "constrained", k, mesh, dS, mu, rho, dt, "cranknicolson", radial);
        st = initial_state(pde, mesh, k);
        E(i, :) = fluid_errors(mesh, k, st.fluid, dS, true);
    end
    n_fail = n_fail + report(sprintf("transfer (disk) k = %d", k), ns, E, [2 2 2]);
end

%% 2. polynomial reproduction
for k = 1:3
    if k == 1
        Up = "(1 + t + t^2)*(x^2 + x*y + y^2)";  Pp = "(1 + 2*t)*(x - 2*y)";
    else
        Up = sprintf("(1 + t + t^2)*(x^%d + x^2*y^%d + y^%d)", k + 1, k - 1, k + 1);
        Pp = "(1 + 2*t)*(x^2 - x*y + y)";
    end
    dpS = IncompressibleElasticityProblem.get_manufactured_solution(Up, Pp, mu, rho, 0);
    dpN = NavierStokesProblem.get_manufactured_solution(Up, Pp, mu, [], rho, 0);
    % the Navier-Stokes source has degree 2k-1: past k = 2 the element
    % rule of test_elements (exact to degree 2k+1) no longer integrates it
    arms = ["stokes", "ns"];
    if k == 3, arms = "stokes"; end
    for arm = arms
        for closure = ["constrained", "pinned"]
            d = ternary(arm == "stokes", dpS, dpN);
            [pde, mesh] = fluid_problem(arm, closure, k, StructuredSimplicialGrid(4, 4).add_stabilization(1), ...
                d, mu, rho, dt, "cranknicolson");
            e = fluid_errors(mesh, k, getfield(initial_state(pde, mesh, k), "fluid"), d, false);
            n_fail = n_fail + check_poly(sprintf("%s %s k = %d", arm, closure, k), e);
        end
    end
    if k < 3
        mesh = immerse_mesh(StructuredSimplicialGrid(8, 8, -1, 1, -1, 1), phi).add_stabilization(1);
        [pde, mesh] = fluid_problem("stokes", "constrained", k, mesh, dpS, mu, rho, dt, "cranknicolson", radial);
        e = fluid_errors(mesh, k, getfield(initial_state(pde, mesh, k), "fluid"), dpS, true);
        n_fail = n_fail + check_poly(sprintf("transfer (disk) k = %d", k), e);
    end
end

%% 3. theta = 1/2 from initial_state vs from the exact Pi u_t(0)
data_at = arrayfun(@(n) IncompressibleElasticityProblem.get_manufactured_solution(U, P, mu, rho, n*dt), 1:4);
for k = 1:2
    ns = [4 8 16];
    D = zeros(numel(ns), 2);
    for i = 1:numel(ns)
        [pde, mesh] = fluid_problem("stokes", "constrained", k, ...
            StructuredSimplicialGrid(ns(i), ns(i)).add_stabilization(1), dS, mu, rho, dt, "cranknicolson");
        s_c = initial_state(pde, mesh, k);
        s_e = s_c;
        s_e.fluid.acceleration = L2proj2D(mesh, dS.u_t, k);
        for n = 1:4
            pde.Data.source_term = data_at(n).f;
            pde.Data.datum_dirichlet = data_at(n).u;
            [s_c, u_c] = step(pde, mesh, k, s_c);
            [s_e, u_e] = step(pde, mesh, k, s_e);
            if n == 1 || n == 4
                D(i, 1 + (n == 4)) = norm(u_c - u_e, "fro")/norm(u_e, "fro");
            end
        end
    end
    fprintf("theta = 1/2, k = %d: |u^n(initial_state) - u^n(exact Pi u_t)| / |u^n|\n", k);
    fprintf("   N    n = 1      EOC     n = 4      EOC\n");
    eoc = compute_eoc(D, 1./ns);
    for i = 1:numel(ns)
        r = NaN(1, 2);
        if i > 1, r = eoc(i-1, :); end
        fprintf("  %3d   %.3e  %5.2f   %.3e  %5.2f\n", ns(i), D(i, 1), r(1), D(i, 2), r(2));
    end
    if any(diff(D) >= 0, "all")
        fprintf("theta = 1/2 k = %d: the difference does not shrink with h FAILED\n", k); n_fail = n_fail + 1;
    end
end

%% 4. laziness
k = 2;
mesh = StructuredSimplicialGrid(8, 8).add_stabilization(1);
pde = fluid_problem("stokes", "constrained", k, mesh, dS, mu, rho, dt, "euler");
profile clear; profile on;
st_be = initial_state(pde, mesh, k);
profile off;
info = profile("info");
built = any(contains({info.FunctionTable.FunctionName}, "AccelerationMatrices"));
lazy_ok = isempty(st_be.fluid.acceleration) && isempty(st_be.fluid.pressure) ...
    && isempty(st_be.fluid.stress) && ~built;
st_wp = initial_state(pde, mesh, k, with_pressure = true);
pde.Data.time_scheme = "cranknicolson";
st_cn = initial_state(pde, mesh, k);
forced_ok = ~isempty(st_wp.fluid.pressure) && isequal(st_wp, st_cn);
fprintf("laziness: theta = 1 -> [] and no AccelerationMatrices call: %d; with_pressure = true " + ...
    "== theta = 1/2 solve: %d\n", lazy_ok, forced_ok);
if ~(lazy_ok && forced_ok)
    fprintf("laziness FAILED\n"); n_fail = n_fail + 1;
end

%% 5. round trip, bitwise; the advanced pressure is -tr(sigma)/2
for scheme = ["euler", "cranknicolson"]
    [pde, mesh] = fluid_problem("stokes", "pinned", k, StructuredSimplicialGrid(8, 8).add_stabilization(1), ...
        dS, mu, rho, dt, scheme);
    s0 = initial_state(pde, mesh, k, with_pressure = true);
    pde.Data.source_term = data_at(1).f;
    pde.Data.datum_dirichlet = data_at(1).u;
    pde.Data.datum_neumann = {data_at(1).sigma{1,1}; data_at(1).sigma{2,1}};
    s1 = step(pde, mesh, k, s0);
    % the same step from a hand-built state
    hand = struct("time", 0, "step", 0, "solid", struct.empty, "fluid", ...
        struct("velocity", L2proj2D(mesh, dS.u, k), "acceleration", s0.fluid.acceleration), ...
        "diagnostics", struct.empty);
    h1 = step(pde, mesh, k, hand);
    same = isequal(s1.fluid.velocity, h1.fluid.velocity) ...
        && isequal(s1.fluid.acceleration, h1.fluid.acceleration) ...
        && isequal(s1.fluid.stress, h1.fluid.stress) ...
        && isequal(s1.fluid.pressure, h1.fluid.pressure) ...
        && s1.time == dt && s1.step == 1;
    fprintf("round trip %s: state == hand-built, bitwise: %d\n", scheme, same);
    if ~same
        fprintf("round trip %s FAILED\n", scheme); n_fail = n_fail + 1;
    end
end
% polynomial data: one Crank-Nicolson step is exact, so -tr(sigma^1)/2 is p(dt)
dp0 = IncompressibleElasticityProblem.get_manufactured_solution( ...
    "(1 + t + t^2)*(x^3 + x^2*y + y^3)", "(1 + 2*t)*(x^2 - x*y + y)", mu, rho, 0);
dp1 = IncompressibleElasticityProblem.get_manufactured_solution( ...
    "(1 + t + t^2)*(x^3 + x^2*y + y^3)", "(1 + 2*t)*(x^2 - x*y + y)", mu, rho, dt);
[pde, mesh] = fluid_problem("stokes", "constrained", k, StructuredSimplicialGrid(4, 4).add_stabilization(1), ...
    dp0, mu, rho, dt, "cranknicolson");
s0 = initial_state(pde, mesh, k);
pde.Data.source_term = dp1.f;
pde.Data.datum_dirichlet = dp1.u;
s1 = step(pde, mesh, k, s0);
p_ex = L2proj2D(mesh, dp1.p, k);
e_p = norm(s1.fluid.pressure - p_ex, "fro")/norm(p_ex, "fro");
fprintf("advanced pressure -tr(sigma^1)/2 vs Pi p(dt), polynomial data: %.2e\n", e_p);
if e_p > 1e-8
    fprintf("advanced pressure FAILED\n"); n_fail = n_fail + 1;
end

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_initial_pressure_fluid:failed", "%d checks failed", n_fail);
end

%% ------------------------------------------------------------------------
function [pde, mesh] = fluid_problem(arm, closure, k, mesh, d, mu, rho, dt, scheme, mapping)
% the manufactured fluid at t = 0; pinned: traction face x = 1 with
% q_N(0) = sigma(0) n, n = (1, 0)
arguments
    arm, closure, k, mesh, d, mu, rho, dt, scheme
    mapping = function_handle.empty
end
gN = {@(x, y) 0*x; @(x, y) 0*x};
if closure == "pinned"
    mesh = apply_tag(mesh, "Neumann", @(X) isapprox(X(1, :), 1));
    gN = {d.sigma{1,1}; d.sigma{2,1}};
end
args = {"viscosity", mu, "density", rho, "time_step", dt, "time_scheme", scheme, ...
    "pressure_mean", closure, "initial_velocity", d.u, "source_term", d.f, ...
    "datum_dirichlet", d.u, "datum_neumann", gN, "datum_dirichlet_rate", d.u_t, ...
    "transfer_mapping", mapping};
if arm == "stokes"
    pde = IncompressibleElasticityProblem(k, args{:});
else
    pde = NavierStokesProblem(k, args{:}, "convective_velocity", d.beta);
end
end

function [state, u] = step(pde, mesh, k, state)
pde.Data.state = state;
sol = solve(pde, mesh, k);
state = get_kinematic_state(pde, sol, state);
u = state.fluid.velocity;
end

function e = fluid_errors(mesh, k, fl, d, shift_mean)
% relative L2 errors [u_t0, p0, sigma0] against the P_k projections
% (sigma componentwise, in PEERS dofs with zero bubbles); shift_mean
% compares p up to its mean on the mesh domain
d2 = nchoosek(k + 2, 2);
N = size(mesh.elements, 2);
M = mass_matrix(mesh, k);
MP = PEERS_matrix(mesh, k);
nrm = @(Mk, v) sqrt(sum(pagemtimes(reshape(v, 1, [], N), pagemtimes(Mk, reshape(v, [], 1, N))), "all"));
vnrm = @(v) sqrt(nrm(M, v(1:d2, :))^2 + nrm(M, v(d2+1:end, :))^2);
a_ex = L2proj2D(mesh, d.u_t, k);
p_ex = L2proj2D(mesh, d.p, k);
s_ex = [L2proj2D(mesh, {d.sigma{1,1}; d.sigma{2,1}; d.sigma{1,2}; d.sigma{2,2}}, k);
        zeros(size(fl.stress, 1) - 4*d2, N)];
ep = fl.pressure - p_ex;
es = fl.stress - s_ex;
if shift_mean
    one = L2proj2D(mesh, @(x, y) 1 + 0*x, k);
    mean_of = @(v) sum(pagemtimes(reshape(one, 1, [], N), pagemtimes(M, reshape(v, [], 1, N))), "all") ...
        / sum(mesh.volume);
    c = mean_of(ep);
    ep = ep - c*one;
    es(1:d2, :) = es(1:d2, :) + c*one;
    es(3*d2 + (1:d2), :) = es(3*d2 + (1:d2), :) + c*one;
    p_ex = p_ex - mean_of(p_ex)*one;
end
e = [vnrm(fl.acceleration - a_ex)/vnrm(a_ex), nrm(M, ep)/nrm(M, p_ex), nrm(MP, es)/nrm(MP, s_ex)];
end

function n_fail = report(name, ns, E, observed)
fprintf("%s\n   N    |u_t0-Pi u_t|  EOC   |p0-Pi p|  EOC   |s0-Pi s|  EOC\n", name);
eoc = compute_eoc(E, 1./ns);
for i = 1:numel(ns)
    r = NaN(1, 3);
    if i > 1, r = eoc(i-1, :); end
    fprintf("  %3d   %.3e  %5.2f   %.3e  %5.2f   %.3e  %5.2f\n", ns(i), E(i, 1), r(1), E(i, 2), r(2), E(i, 3), r(3));
end
n_fail = double(any(r < observed - 0.3));
if n_fail
    fprintf("%s: EOC below the measured order - 0.3 (%s) FAILED\n", name, mat2str(observed));
end
end

function n_fail = check_poly(name, e)
fprintf("polynomial %-24s |u_t0| %.2e  |p0| %.2e  |sigma0| %.2e\n", name, e);
n_fail = double(max(e) > 1e-8);
if n_fail, fprintf("polynomial %s FAILED\n", name); end
end

function v = ternary(c, a, b)
if c, v = a; else, v = b; end
end
