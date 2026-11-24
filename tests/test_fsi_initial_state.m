% test_fsi_initial_state.m
% FluidStructureCoupling.initial_state: the t = 0 state of the coupled
% problem (solid state -> transported solid acceleration on Sigma -> fluid
% acceleration problem with the kinematic datum a_F = a_S + d/dt T2).
%  1. Rates on the flat manufactured problem
%       u_F = curl(sin(pi x) sin(pi y) cos(t + x y)),
%       p   = sin(pi x) exp(y - 1/2) (1 + t),
%       u_S = (x sin(y) cos(t), y cos(x) (1 + t)),
%     rho_F = rho_S = 1, mu = 0.4, E = 1, nu = 0.25, newmark solid and
%     Crank-Nicolson fluid, against the P_k projections at t = 0 (relative
%     L2): the solid a0, a_S on Sigma (the transported datum at the Sigma
%     quadrature points, discrete L2(Sigma)), the fluid u_t0, p0
%     (p = -tr(sigma_F)/2) and sigma_F0, and the relative traction
%     residual. Families: fitted (fitted_split, gap 0; pinned closure
%     k = 1, 2, 3 and constrained k = 2) and unfitted (Sigma = {y = 1/2}
%     carved out of the "left" N x N grid, N = 9, 17, 33, cut elements
%     discarded: both meshes retreat by h/2; pinned, k = 1, 2). The
%     newmark solid starts from the HDG elliptic projection of its data,
%     so its a0 is the projection P_k u_tt(0) itself: asserted at
%     round-off (1e-10 relative, "r.o." below). The other columns:
%     EOC >= asserted - 0.3 on the last refinement:
%                             a0    a_S(Sigma)  u_t0  p0   sigma0  residual
%       fitted   k = 1, 2, 3  r.o.  k+1         k+1   k+1  k+1     k+1
%       unfitted k = 1        r.o.  2           1     2    2       1.8
%       unfitted k = 2        r.o.  3           2     1.5  1.5     1.5
%     (unfitted, observed at 1/h = 17, 33: k = 1 u_t0 1.94, 1.96,
%     p0 1.88, 1.95, sigma0 1.86, 1.94, residual 2.05, 2.01; k = 2 u_t0
%     2.80, 2.76, p0 3.68, 3.47, sigma0 3.65, 3.45, residual 3.12, 3.16;
%     the unfitted bounds of those columns stay below them because the
%     finest mesh is still coarse).
%  2. Polynomial data (test_fsi_transient_polynomial: degree <= k,
%     g(t) = 1 + t + t^2/2, u_F = B (1 + t)) at gap 0 and gap 0.02, k = 2:
%     a0, its skeleton trace, u_t0, p0 reproduce the projections and the
%     traction residual vanishes (1e-8 relative), and one coupled
%     Newmark / Crank-Nicolson step from the state is exact (1e-10).
%  3. The incompatibility flag: compatible data raise no warning; the
%     solid Young modulus doubled after load_problems (a traction mismatch
%     of the order of the traction) and a rest start under a solid load
%     raise FluidStructureCoupling:incompatibleInitialData.
%  4. Laziness: a backward-Euler fluid without with_pressure solves no
%     acceleration problem (no AccelerationMatrices call; acceleration,
%     pressure, stress and the traction residual are []); with_pressure =
%     true solves it.
%  5. A onestep solid with a Crank-Nicolson fluid is rejected
%     (FluidStructureCoupling:onestepAcceleration).
clc; clear;
n_fail = 0;
E = 1.0; nu = 0.25; mu = 0.4; rho = 1.0; dt = 1e-2;
UU = "sin(pi*x)*sin(pi*y)*cos(t + x*y)";
PP = "sin(pi*x)*exp(y - 0.5)*(1 + t)";
US = ["x*sin(y)*cos(t)"; "y*cos(x)*(1 + t)"];
smooth = {UU, PP, US};

%% 1. rates
flat = FlatCurveHorizontal();
arms = struct("name", {}, "k", {}, "closure", {}, "unfitted", {}, "hs", {}, "observed", {});
fitted_obs = @(k) [NaN, k + 1, k + 1, k + 1, k + 1, k + 1];   % NaN: round-off
for k = 1:3
    arms(end + 1) = struct("name", "fitted", "k", k, "closure", "pinned", ...
        "unfitted", false, "hs", 1./[4 8 16 32], "observed", fitted_obs(k)); %#ok<SAGROW>
end
arms(end + 1) = struct("name", "fitted", "k", 2, "closure", "constrained", ...
    "unfitted", false, "hs", 1./[4 8 16 32], "observed", fitted_obs(2));
arms(end + 1) = struct("name", "unfitted", "k", 1, "closure", "pinned", ...
    "unfitted", true, "hs", 1./[9 17 33], "observed", [NaN, 2, 1, 2, 2, 1.8]);
arms(end + 1) = struct("name", "unfitted", "k", 2, "closure", "pinned", ...
    "unfitted", true, "hs", 1./[9 17 33], "observed", [NaN, 3, 2, 1.5, 1.5, 1.5]);
for arm = arms
    k = arm.k;
    Err = zeros(numel(arm.hs), 6);
    for i = 1:numel(arm.hs)
        if arm.unfitted
            N = round(1/arm.hs(i));
            [fm, sm] = StructuredSimplicialGrid(N, N, "left").split(flat, cut_elements = "discard");
        else
            [fm, sm] = fitted_split(flat, arm.hs(i), 0);
        end
        Err(i, :) = t0_errors(k, flat, fm, sm, smooth, rho, mu, E, nu, dt, arm.closure);
    end
    n_fail = n_fail + report(sprintf("%s %s k = %d", arm.name, arm.closure, k), arm.hs, Err, arm.observed);
end

%% 2. polynomial data, gap 0 and gap 0.02
k = 2;
poly = {"x*y*(1 + t)", "(x + y)*(1 + t)", ["x*(1 + t + t^2/2)"; "y*(1 + t + t^2/2)"]};
for gap = [0, 0.02]
    [fm, sm] = fitted_split(flat, 1/8, gap);
    [fp, sp, src, d] = coupled_problems(k, flat, 0.0, poly, rho, mu, E, nu, dt, "pinned");
    st = FluidStructureCoupling.initial_state(fp, sp, fm, sm, flat, src{:}, ...
        source_dirichlet_rate = d.T_2_t);
    e = [rel(st.solid.acceleration, L2proj2D(sm, d.u_tt_S, k)), ...
         rel(st.solid.acceleration_skeleton, L2projSkeleton1d(sm, d.u_tt_S, k)), ...
         rel(st.fluid.acceleration, L2proj2D(fm, d.u_t_F, k)), ...
         rel(st.fluid.pressure, L2proj2D(fm, {@(x, y) -(d.sigma_F{1,1}(x, y) + d.sigma_F{2,2}(x, y))/2}, k)), ...
         st.diagnostics.traction_residual_relative];
    [fp1, sp1, src1, d1] = coupled_problems(k, flat, dt, poly, rho, mu, E, nu, dt, "pinned");
    [fs, ss] = FluidStructureCoupling.solve(fp1, sp1, fm, sm, flat, src1{:}, "newmark", state = st);
    es = ss.compute_l2_error({d1.sigma_S, d1.u_S, d1.gamma_S}, {d1.u_S});
    ef = fs.compute_l2_error({d1.sigma_F, d1.u_F, d1.gamma_F}, {d1.u_F});
    step_error = max([es(:); ef(:)]);
    fprintf("polynomial gap %-4g |a0| %.2e  |ah0| %.2e  |u_t0| %.2e  |p0| %.2e  residual %.2e  " + ...
        "one step %.2e\n", gap, e, step_error);
    if max(e) > 1e-8 || step_error > 1e-10
        fprintf("polynomial gap %g FAILED\n", gap); n_fail = n_fail + 1;
    end
end

%% 3. the incompatibility flag
k = 1;
[fm, sm] = fitted_split(flat, 1/4, 0);
for variant = ["compatible", "Young modulus doubled", "rest under a load"]
    [fp, sp, src, d] = coupled_problems(k, flat, 0.0, smooth, rho, mu, E, nu, dt, "pinned");
    rate = d.T_2_t;
    switch variant
        case "Young modulus doubled"
            sp.Data.young_modulus = 2*E;
        case "rest under a load"
            zero = {@(x, y) 0*x; @(x, y) 0*x};
            [sp.Data.initial_displacement, sp.Data.initial_velocity, fp.Data.initial_velocity, ...
                fp.Data.source_term, fp.Data.datum_dirichlet, fp.Data.datum_dirichlet_rate, ...
                sp.Data.datum_dirichlet, rate] = deal(zero);
            sp.Data.source_term = {@(x, y) 0*x; @(x, y) -2 + 0*x};
            src = {zero, zero};
    end
    lastwarn("");
    st = FluidStructureCoupling.initial_state(fp, sp, fm, sm, flat, src{:}, ...
        source_dirichlet_rate = rate);
    [~, id] = lastwarn;
    flagged = id == "FluidStructureCoupling:incompatibleInitialData";
    fprintf("flag, %s: relative traction residual %.3e, warning %d\n", ...
        variant, st.diagnostics.traction_residual_relative, flagged);
    if flagged ~= (variant ~= "compatible")
        fprintf("flag FAILED\n"); n_fail = n_fail + 1;
    end
end

%% 4. laziness
k = 2;
[fm, sm] = fitted_split(flat, 1/8, 0);
[fp, sp, src, d] = coupled_problems(k, flat, 0.0, smooth, rho, mu, E, nu, dt, "pinned");
fp.Data.time_scheme = "euler";
profile clear; profile on;
st = FluidStructureCoupling.initial_state(fp, sp, fm, sm, flat, src{:}, source_dirichlet_rate = d.T_2_t);
profile off;
info = profile("info");
built = any(contains({info.FunctionTable.FunctionName}, "AccelerationMatrices"));
lazy_ok = ~built && isempty(st.fluid.acceleration) && isempty(st.fluid.pressure) ...
    && isempty(st.fluid.stress) && isempty(st.diagnostics.traction_residual) ...
    && ~isempty(st.fluid.velocity) && ~isempty(st.solid.acceleration);
st = FluidStructureCoupling.initial_state(fp, sp, fm, sm, flat, src{:}, ...
    source_dirichlet_rate = d.T_2_t, with_pressure = true);
forced_ok = ~isempty(st.fluid.acceleration) && ~isempty(st.fluid.pressure) ...
    && ~isempty(st.diagnostics.traction_residual);
fprintf("laziness: euler -> no acceleration problem: %d; with_pressure = true solves it: %d\n", ...
    lazy_ok, forced_ok);
if ~(lazy_ok && forced_ok)
    fprintf("laziness FAILED\n"); n_fail = n_fail + 1;
end

%% 5. onestep solid + Crank-Nicolson fluid
sp.Data.time_scheme = "onestep";
fp.Data.time_scheme = "cranknicolson";
try
    FluidStructureCoupling.initial_state(fp, sp, fm, sm, flat, src{:});
    id = "";
catch err
    id = err.identifier;
end
fprintf("onestep solid + Crank-Nicolson fluid raises: %s\n", id);
if id ~= "FluidStructureCoupling:onestepAcceleration"
    fprintf("onestep rejection FAILED\n"); n_fail = n_fail + 1;
end

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_fsi_initial_state:failed", "%d checks failed", n_fail);
end

%% ------------------------------------------------------------------------
function [fp, sp, src, d] = coupled_problems(k, I, t, data, rho, mu, E, nu, dt, closure)
% load_problems at time t, newmark solid and Crank-Nicolson fluid
[fp, sp, src, d] = FluidStructureCoupling.load_problems(k, I, t, dt, data{:}, ...
    rho, mu, rho, E, nu, 0, pressure_mean = closure);
sp.Data.time_scheme = "newmark";
fp.Data.time_scheme = "cranknicolson";
end

function e = t0_errors(k, I, fm, sm, data, rho, mu, E, nu, dt, closure)
% relative errors [a0, a_S on Sigma, u_t0, p0, sigma_F0] and the relative
% traction residual
[fp, sp, src, d] = coupled_problems(k, I, 0.0, data, rho, mu, E, nu, dt, closure);
st = FluidStructureCoupling.initial_state(fp, sp, fm, sm, I, src{:}, source_dirichlet_rate = d.T_2_t);
d2 = nchoosek(k + 2, 2);
N = size(fm.elements, 2);
M = mass_matrix(fm, k);
MP = PEERS_matrix(fm, k);
nrm = @(Mk, v) sqrt(sum(pagemtimes(reshape(v, 1, [], N), pagemtimes(Mk, reshape(v, [], 1, N))), "all"));
vnrm = @(v) sqrt(nrm(M, v(1:d2, :))^2 + nrm(M, v(d2+1:end, :))^2);
a_ex = L2proj2D(fm, d.u_t_F, k);
p_ex = L2proj2D(fm, {@(x, y) -(d.sigma_F{1,1}(x, y) + d.sigma_F{2,2}(x, y))/2}, k);
s_ex = [L2proj2D(fm, {d.sigma_F{1,1}; d.sigma_F{2,1}; d.sigma_F{1,2}; d.sigma_F{2,2}}, k);
        zeros(size(st.fluid.stress, 1) - 4*d2, N)];
e_fluid = [vnrm(st.fluid.acceleration - a_ex)/vnrm(a_ex), ...
           nrm(M, st.fluid.pressure - p_ex)/nrm(M, p_ex), ...
           nrm(MP, st.fluid.stress - s_ex)/nrm(MP, s_ex)];
% the solid acceleration, in the volume and transported to Sigma
Ms = mass_matrix(sm, k);
Ns = size(sm.elements, 2);
snrm = @(v) sqrt(sum(pagemtimes(reshape(v(1:d2, :), 1, [], Ns), pagemtimes(Ms, reshape(v(1:d2, :), [], 1, Ns))), "all") ...
               + sum(pagemtimes(reshape(v(d2+1:end, :), 1, [], Ns), pagemtimes(Ms, reshape(v(d2+1:end, :), [], 1, Ns))), "all"));
a_s_ex = L2proj2D(sm, d.u_tt_S, k);
e_a = snrm(st.solid.acceleration - a_s_ex)/snrm(a_s_ex);
[~, p2] = I.get_partitions(fm, sm, quadrature_patch = k + 1, quadrature_path = k + 2);
ext = trace_extension_coefficients(sp, sm, p2, st.solid.stress_acceleration, ...
    st.solid.rotation_acceleration, Order = k);
a_sigma = interface_datum(I, fm, p2, st.solid.acceleration_skeleton(:, p2.faces(:)) + ext);
X = p2.quadrature.interface;
x = squeeze(X(1, :, :));  y = squeeze(X(2, :, :));
[~, w] = lgwt(k + 1, 0, 1);
W = w(:) .* reshape(p2.area, 1, []);
diff2 = 0;  ref2 = 0;
for c = 1:2
    diff2 = diff2 + sum(W .* (a_sigma{c}(x, y) - d.u_tt_S{c}(x, y)).^2, "all");
    ref2 = ref2 + sum(W .* d.u_tt_S{c}(x, y).^2, "all");
end
e = [e_a, sqrt(diff2/ref2), e_fluid, st.diagnostics.traction_residual_relative];
end

function r = rel(a, b)
r = norm(a - b, "fro")/norm(b, "fro");
end

function n_fail = report(name, hs, Err, observed)
fprintf("%s\n    1/h   |a0|      EOC   |a_S(Sigma)| EOC  |u_t0|     EOC   |p0|      EOC   " + ...
    "|sigma0|  EOC   residual  EOC\n", name);
eoc = compute_eoc(Err, hs);
for i = 1:numel(hs)
    r = NaN(1, 6);
    if i > 1, r = eoc(i-1, :); end
    fprintf("  %4d  %s\n", round(1/hs(i)), sprintf("%.3e %5.2f   ", [Err(i, :); r]));
end
round_off = isnan(observed);
n_fail = double(any(r(~round_off) < observed(~round_off) - 0.3) ...
    || any(max(Err(:, round_off), [], 1) > 1e-10));
if n_fail
    fprintf("%s: EOC below the asserted order - 0.3, or a round-off column above 1e-10 " + ...
        "(%s, NaN = round-off) FAILED\n", name, mat2str(observed));
end
end
