function out = eulerian_hyperelastic_fsi_march(k, N, amplitude, dt, n_steps, uu_fluid, p_fluid, options)
%EULERIAN_HYPERELASTIC_FSI_MARCH  FSI march with the fully Eulerian hyperelastic solid, both sides carved.
%
%   out = eulerian_hyperelastic_fsi_march(k, N, A, dt, n_steps, uu_fluid, p_fluid, ...
%       rotation = 0, law = "linear", w = [...])
%
%   Sigma(t) is the circle of the given radius and centre c(t) = center +
%   A sin(2 pi t) e_x; at every t_n both meshes are the two sides of
%   StructuredSimplicialGrid(N, N, "crossed").split by Sigma(t_n), cut
%   elements discarded, interface vertices snapped onto Sigma(t_n).
%   geometry = "flat": Sigma(t) = {y = y0 + A sin(2 pi t)}, solid below
%   (translating_flat_problems' fluid data); the solid's walls (its
%   boundary faces off Sigma) carry u_hat and v_hat = the face L2
%   projections of the exact u and v (their trace rows replaced).
%   Fluid: Stokes, PEERS stress, velocity and pressure, the data of
%   translating_circle_problems (uu_fluid, p_fluid); its two velocity
%   levels are carried across carves by fluid_transfer = "remap"
%   (velocity_stress_fsi_march's rule: an element that switches on or
%   changes shape takes the L2 projection of one old element's
%   polynomial) or "projection" (the solid's zero-shift two-mesh L2
%   projection with the extension of extension_degree).
%   Solid: the fully Eulerian scheme of eulerian_hyperelastic_solid_march
%   (unknowns sigma, u, gamma, v; u_hat, v_hat; law, rotation, w as there,
%   fields "standard"; Newton on the transport and the stress), its levels
%   carried by the zero-shift two-mesh projection, the displacement
%   skeleton closed on Sigma by the characteristic rule with no
%   displacement data: displacement_data = "kinematic" (default), g_u the
%   BDF2 step of the kinematic equation at the face points, Z_u + (2 dt/3)
%   (v - (grad u) v), from the iterate (in the Newton Jacobian);
%   "extrapolated", g_u the trace of 2 X~^n_u - X~^n-1_u (an error of
%   order dt^2 |d_tt u| in u_hat on inflow faces, amplified by about k^2/h
%   in L_h: the solid errors stop converging in h for k >= 2).
%   One BDF2 step on both sides from the first step (the two initial
%   levels: the exact fields at t = 0 and t = -dt, the fluid's projected
%   on the fluid carve of t = 0, the solid's on the solid carves of t = 0
%   and -dt); the fluid is the theta = 2/3 step from (4 u^n - u^n-1)/3.
%   Each Newton iteration is one coupled solve by
%   FluidStructureCoupling.interface_problem ("onestep", time step dt,
%   swap_transmission) with the solid given as Data.eulerian_blocks: the
%   solid's element blocks at the iterate in the unknown w = dt v (so that
%   the kinematic row <dt u_F - w_S, mu> = dt <T2, mu> carries no history),
%   the solid mesh tau = tau_v/dt (the T1 row's solid pair tau_v (v -
%   v_hat)). Transmission data from the exact fields on each side: T1 =
%   sigma_F n_F + sigma_S n_S, T2 = u_F - v_S (the solid's Cauchy stress
%   and velocity of eulerian_hyperelastic_solid_march's exact fields).
%   interface = "solid" (circle only): Sigma is not prescribed but taken
%   from the solid, Sigma = {x : x - u(x) in Sigma_0}, Sigma_0 the circle
%   about center (eulerian_interface_curve): at t = 0 and -dt from the
%   exact u, and for the carve of t_n+1 from the solid's u^n on the carve
%   of t_n (lagged one step; interface_predictor = "velocity" moves the
%   material points by dt v^n); T1 then uses the normal of that curve.
%   out.interface_samples (n_steps x 4): the samples of Sigma evaluated
%   outside the solid mesh, the largest |x - X| on Sigma, and the largest
%   distance (at equal reference parameter) of the carved Sigma from the
%   exact Sigma of t_n (the reconstruction error) and of t_n+1 (with the
%   lag).
%   convection = true: Navier-Stokes fluid; the fluid source gains rho_f
%   (u_F.grad) u_F and every coupled solve convects with beta = u*, the
%   divergence-free postprocessed velocity (divergence_free_velocity) of
%   the previous iterate (at the first iterate of a step, the last u*
%   remapped to the new carve at degree k + 1; at the first step, none:
%   Stokes), so that each Newton iteration is also a Picard step of the
%   convection; the loop stops only after a convected solve.
%   out.divergence_star (n_steps): the largest max_K ||div u*|| of a step.
%   operators = true: no march; out holds handles to declarative_blocks,
%   wall_rows, remap_fluid_state and vector_error.
%   The loop stops when the relative change of the solid unknowns (v, u,
%   v_hat, u_hat) and the fluid velocity is at most tolerance.
%   swap_transmission (default true) as in interface_problem; solid_tau
%   ([] = tau_v) replaces tau_v in the T1 row's solid pair.
%
%   out: step_errors (n_steps x 4: [u_F, v_S, u_S, sigma_S] L2 errors),
%   step_norms (the exact norms), iterations, converged, deviation_energy
%   (the solid's, as in eulerian_hyperelastic_solid_march); the final v,
%   u, sigma (solid element dofs) and fluid_velocity, the final meshes;
%   seconds.
%   initial_perturbation = [seed, amplitude]: both initial solid levels
%   get a random (dv, du) of linearized energy amplitude^2
%   (out.perturbation_energy, its energy); record_state = true:
%   out.record(n) = struct(step, level, lam, fluid_velocity) after each
%   step; reference = such a record of a run on the same carves:
%   out.difference (n_steps x 2) = [the solid's linearized energy (as in
%   eulerian_hyperelastic_solid_march), rho_f/2 |u_F|^2] of this run's
%   state minus the reference's.
%
%See also: eulerian_hyperelastic_solid_march, velocity_stress_fsi_march,
%   FluidStructureCoupling.interface_problem
arguments
    k (1,1) double
    N (1,1) double
    amplitude (1,1) double
    dt (1,1) double
    n_steps (1,1) double
    uu_fluid (1,1) string
    p_fluid (1,1) string
    options.rotation (1,1) double = 0
    options.law (1,1) string = "linear"
    options.w (2,1) string = ["0.02*sin(5*X + 3*Y)*cos(3*t)"; "0.02*cos(4*X - 5*Y)*cos(2*t)"]
    options.extension_degree double = []
    options.geometry (1,1) string {mustBeMember(options.geometry, ["circle", "flat"])} = "circle"
    options.y0 (1,1) double = 0.45
    options.center (2,1) double = [0.5; 0.5]
    options.radius (1,1) double = 0.2
    options.rho_f (1,1) double = 1
    options.mu (1,1) double = 1
    options.rho_s (1,1) double = 1
    options.E (1,1) double = 10
    options.nu (1,1) double = 0.3
    options.tolerance (1,1) double = 1e-10
    options.max_iterations (1,1) double = 20
    options.swap_transmission (1,1) logical = true
    options.solid_tau double = []
    options.fluid_transfer (1,1) string {mustBeMember(options.fluid_transfer, ["remap", "projection"])} = "remap"
    options.initial_perturbation double = []
    options.record_state (1,1) logical = false
    options.reference struct = struct.empty
    options.displacement_data (1,1) string {mustBeMember(options.displacement_data, ...
        ["kinematic", "extrapolated"])} = "kinematic"
    options.interface (1,1) string {mustBeMember(options.interface, ["prescribed", "solid"])} = "prescribed"
    options.interface_predictor (1,1) string {mustBeMember(options.interface_predictor, ...
        ["none", "velocity"])} = "none"
    options.convection (1,1) logical = false
    options.operators (1,1) logical = false
end
if options.operators
    % the pieces another coupled march builds on (the positionals unused)
    out = struct("declarative_blocks", @declarative_blocks, "wall_rows", @wall_rows, ...
        "remap_fluid_state", @remap_fluid_state, "vector_error", @vector_error);
    return
end
wall = tic;
ops = eulerian_hyperelastic_solid_march(k, N, amplitude, dt, 0, rotation = options.rotation, ...
    law = options.law, w = options.w, displacement_boundary = "characteristic", ...
    extension_degree = options.extension_degree, center = options.center, radius = options.radius, ...
    rho = options.rho_s, E = options.E, nu = options.nu, operators = true);
prm = ops.prm;
ex = ops.exact;
solid_tau = prm.tau_v;
if ~isempty(options.solid_tau), solid_tau = options.solid_tau; end
d1 = k + 1;
d2 = nchoosek(k + 2, 2);
background = StructuredSimplicialGrid(N, N, "crossed");
physics = {options.rho_f, options.mu, options.rho_s, options.E, options.nu};
flat = options.geometry == "flat";
if flat
    curve_at = @(t) FlatCurveHorizontal(options.y0 + amplitude*sin(2*pi*t));
    normal_S = @(x, y, curve) {0*x, 1 + 0*x};
    data = @(t) fluid_data_flat(k, options.y0, amplitude, t, dt, uu_fluid, p_fluid, physics);
else
    curve_0 = CircularCurve(options.center, options.radius);
    curve_at = @(t) CircularCurve(options.center + [amplitude*sin(2*pi*t); 0], options.radius);
    normal_S = @(x, y, curve) {(x - curve.center(1))/curve.radius, (y - curve.center(2))/curve.radius};
    data = @(t) fluid_data(k, curve_0, amplitude, t, dt, uu_fluid, p_fluid, physics);
end
from_solid = options.interface == "solid";
if from_solid
    % Sigma(t) = {x : x - u(x) in Sigma_0}, Sigma_0 the circle about the
    % reference centre: from the exact u at t = 0 and -dt, then from the
    % solid's u^n for the carve of t_n+1
    assert(~flat, "eulerian_hyperelastic_fsi_march:interface", "interface = ""solid"" needs the circle.");
    exact_u = @(t) @(x) reshape(ex.u(x(1, :).', x(2, :).', t), [], 2).';
    curve_at = @(t) eulerian_interface_curve(curve_0, [], [], exact_u(t), spacing = 0.25/N);
    normal_S = @(x, y, curve) curve_normal(curve, x, y);
end
interface_samples = zeros(n_steps, 4);   % [outside samples, max |d|, deviation from Sigma(t_n), from Sigma(t_n+1)]
div_star = zeros(n_steps, 1);            % max_K ||div u*|| over the iterates of each step
quad = {"n_quadpoints", k + 2, "n_quadpoints_path", k + 3};
if options.convection
    % Navier-Stokes: the fluid source gains rho_f (u_F.grad) u_F
    conv = convection_source(uu_fluid, options.rho_f);
    data = @(t) with_convection(data(t), conv, t);
end
beta = {};   % the convective velocity u* (P_k+1 dofs) of the last iterate

% initial levels
[fm, sm, keys] = carve_both(background, curve_at(0));
[~, sm_prev] = carve_both(background, curve_at(-dt));
D0 = data(0);
Dm = data(-dt);
UF = project_handles(fm, k, D0.pd.u_F);
UF_prev = project_handles(fm, k, Dm.pd.u_F);
level = ops.project_vector(sm, k, @(x, y) [ex.v(x, y, 0); ex.u(x, y, 0)]);
level_prev = ops.project_vector(sm_prev, k, @(x, y) [ex.v(x, y, -dt); ex.u(x, y, -dt)]);
perturbation_energy = [];
if ~isempty(options.initial_perturbation)
    % a random solid (dv, du) of linearized energy amplitude^2 added to
    % both initial levels (its level -1 copy projected to carve(-dt))
    geo0 = ops.geometry(sm);
    stream = RandStream("mt19937ar", Seed = options.initial_perturbation(1));
    d = randn(stream, size(level));
    none = struct("dir", zeros(0, 1), "lam", zeros(4*geo0.n_dof, 1));
    d = d/sqrt(ops.linearized_energy(geo0, k, prm, d, ops.initial_skeleton(geo0, k, d, none)));
    d = options.initial_perturbation(2)*d;
    level = level + d;
    level_prev = level_prev + ops.project_between(sm, sm_prev, d);
    perturbation_energy = ops.linearized_energy(geo0, k, prm, d, ops.initial_skeleton(geo0, k, d, none));
end
record = struct("step", {}, "level", {}, "lam", {}, "fluid_velocity", {});
difference = nan(n_steps, 2);

step_errors = nan(n_steps, 4);
step_norms = nan(n_steps, 4);
iterations = zeros(n_steps, 1);
converged = false(n_steps, 1);
deviation_energy = nan(n_steps, 1);
for n = 1:n_steps
    t1 = n*dt;
    D1 = data(t1);
    if from_solid
        % lagged: from u^n on the carve of t_n (predictor "velocity": the
        % material points moved by dt v^n)
        predictor = {};
        if options.interface_predictor == "velocity"
            predictor = {"velocity", level(1:2*d2, :), "time_step", dt};
        end
        [curve, ci] = eulerian_interface_curve(curve_0, sm, k, level(2*d2 + 1:end, :), predictor{:});
        % deviation from the exact Sigma of t_n (the reconstruction) and
        % of t_n+1 (reconstruction plus lag), at equal reference parameter
        ts = linspace(curve_0.t_range(1), curve_0.t_range(2), 401);
        P = curve.get_2d_points(ts);
        dev = @(c) max(vecnorm(P - c.get_2d_points(ts)));
        interface_samples(n, :) = [ci.outside, ci.max_step, dev(curve_at(t1 - dt)), dev(curve_at(t1))];
    else
        curve = curve_at(t1);
    end
    [fm_old, sm_old, keys_old] = deal(fm, sm, keys);
    [fm, sm, keys] = carve_both(background, curve);
    if options.fluid_transfer == "remap"
        UF = remap_fluid_state(UF, fm_old, fm, k, keys_old.fluid, keys.fluid);
        UF_prev = remap_fluid_state(UF_prev, fm_old, fm, k, keys_old.fluid, keys.fluid);
    else
        UF = ops.project_between(fm_old, fm, UF);
        UF_prev = ops.project_between(fm_old, fm, UF_prev);
    end
    if ~isempty(beta)
        % the last u* on the new carve (remap at degree k + 1; exactly
        % divergence-free on the elements that kept their shape)
        b = remap_fluid_state([beta{1}; beta{2}], fm_old, fm, k + 1, keys_old.fluid, keys.fluid);
        beta = {b(1:end/2, :); b(end/2 + 1:end, :)};
    end
    Xn = ops.project_between(sm_old, sm, level);
    Xm = ops.project_between(sm_prev, sm, level_prev);
    sm_prev = sm_old;
    geo = ops.geometry(sm);
    if flat
        % the solid's walls (its boundary faces off Sigma): the exact u and
        % v at t_n+1 at the face points of every element
        is_boundary = accumarray(sm.boundaries(:), 1, [size(sm.faces, 2), 1]) == 1;
        is_wall = is_boundary & ~ismember((1:size(sm.faces, 2)).', sm.intfaces(:));
        wall_faces = is_wall(sm.boundaries);
        xyq = sm.get_2d_physical(repmat(geo.q2.', 1, 1, size(sm.elements, 2)), 1:size(sm.elements, 2));
        xq = reshape(xyq(1, :, :), [], size(sm.elements, 2));
        yq = reshape(xyq(2, :, :), [], size(sm.elements, 2));
        g_wall = cat(3, component(ex.u, xq, yq, n*dt, 1), component(ex.u, xq, yq, n*dt, 2), ...
            component(ex.v, xq, yq, n*dt, 1), component(ex.v, xq, yq, n*dt, 2));
    end
    n_el = size(sm.elements, 2);
    % solid data: sources, no Dirichlet dofs, g_u the extrapolated trace
    Z = (4*Xn - Xm)/3;
    level0 = 2*Xn - Xm;
    sdata.f = ops.moments(sm, k, @(x, y) ex.f(x, y, t1));
    sdata.g = ops.moments(sm, k, @(x, y) ex.g(x, y, t1));
    sdata.lam = zeros(4*geo.n_dof, 1);
    sdata.dir = zeros(0, 1);
    Ue = level0(2*d2 + 1:end, :);
    sdata.gu_face = cat(3, geo.Pf*Ue(1:d2, :), geo.Pf*Ue(d2 + 1:end, :));
    if options.displacement_data == "kinematic"
        sdata.kinematic = struct("Zu", Z(2*d2 + 1:end, :), "dt", dt);
    end
    lam0 = ops.initial_skeleton(geo, k, level0, sdata);
    % fluid: BDF2 (theta = 2/3) from (4 u^n - u^n-1)/3
    fp = D1.fp;
    fp.Data.time_scheme = "cranknicolson";
    fp.Data.theta = 2/3;
    fp.Data.time_step = dt;
    fluid_history = (4*UF - UF_prev)/3;
    state = struct("fluid", struct("velocity", fluid_history, "acceleration", zeros(size(fluid_history))), ...
        "solid", struct());
    T1 = cell(2, 1);
    T2 = cell(2, 1);
    for i = 1:2
        T1{i} = @(x, y) traction(ex, D1.pd.sigma_F, normal_S(x, y, curve), x, y, t1, i);
        T2{i} = @(x, y) D1.pd.u_F{i}(x, y) + 0*x - component(ex.v, x, y, t1, i);
    end
    sm_solve = sm.add_stabilization(solid_tau/dt);
    sp = LinearElasticityProblem(k, density = options.rho_s, time_step = dt, ...
        young_modulus = options.E, poisson_ratio = options.nu);
    ok = false;
    UF_it = fluid_history;
    for it = 1:options.max_iterations
        [A1, A2, Af, A3, A4, Ab] = ops.local_operators(geo, k, ops.hp, prm, sdata, Z, level0, lam0);
        if flat
            [A3, A4, Ab] = wall_rows(geo, A3, A4, Ab, wall_faces, g_wall, d1);
        end
        sp.Data.eulerian_blocks = declarative_blocks(A1, A2, Af, A3, A4, Ab, geo.d2peers, d2, d1, dt, n_el);
        convected = options.convection && ~isempty(beta);
        if convected
            fp.Data.convective_velocity = beta;   % u* of the last iterate, never u_h
        end
        coupled = FluidStructureCoupling.interface_problem(fp, sp, fm, sm_solve, curve, T1, T2, "onestep", ...
            "Order", k, "state", state, "patch_weight", "jacobian", quad{:}, ...
            "swap_transmission", options.swap_transmission);
        [fs, ss] = solve(coupled, k);
        S = ss.Unknowns.field(1).dof;
        level_new = [ss.Unknowns.field(2).dof/dt; ss.Unknowns.field(4).dof];
        what = ss.Unknowns.hybrid(1).dof;   % 2 d1 x n_faces, x over y
        uhat = ss.Unknowns.hybrid(2).dof;
        lam_new = [reshape(uhat(1:d1, :), [], 1); reshape(uhat(d1 + 1:end, :), [], 1); ...
            reshape(what(1:d1, :), [], 1)/dt; reshape(what(d1 + 1:end, :), [], 1)/dt];
        UF_new = fs.Unknowns.field(2).dof;
        change = norm([level_new(:) - level0(:); lam_new - lam0; UF_new(:) - UF_it(:)]) ...
            /max(norm([level_new(:); lam_new; UF_new(:)]), realmin);
        [level0, lam0, UF_it] = deal(level_new, lam_new, UF_new);
        if options.convection
            [u_star, divergence_star] = divergence_free_velocity(fs);
            beta = u_star;
            div_star(n) = max(div_star(n), divergence_star);
        end
        if change <= options.tolerance && (convected || ~options.convection)
            ok = true;
            break
        end
    end
    if ~ok
        warning("eulerian_hyperelastic_fsi_march:newton", ...
            "Newton not converged at step %d: relative change %.2e after %d iterations.", n, change, it);
    end
    iterations(n) = it;
    converged(n) = ok;
    [level_prev, level] = deal(level, level0);
    [UF_prev, UF] = deal(UF, UF_it);
    [eF, nF] = vector_error(fm, k, UF, D1.pd.u_F);
    [eS, nS] = ops.errors(geo, k, level, S, ex, t1);
    step_errors(n, :) = [eF, eS];
    step_norms(n, :) = [nF, nS];
    deviation_energy(n) = ops.deviation(geo, k, prm, level, lam0, ex, t1);
    if options.record_state
        record(n) = struct("step", n, "level", level, "lam", lam0, "fluid_velocity", UF);
    end
    if ~isempty(options.reference)
        r = options.reference([options.reference.step] == n);
        dU = UF - r.fluid_velocity;
        zero = {@(x, y) 0*x, @(x, y) 0*x};
        difference(n, :) = [ops.linearized_energy(geo, k, prm, level - r.level, lam0 - r.lam), ...
            options.rho_f/2*vector_error(fm, k, dU, zero)^2];
    end
end
out = struct("step_errors", step_errors, "step_norms", step_norms, "iterations", iterations, ...
    "converged", converged, "deviation_energy", deviation_energy, "v", level(1:2*d2, :), ...
    "u", level(2*d2 + 1:end, :), "sigma", S, "fluid_velocity", UF, "fluid_mesh", fm, "solid_mesh", sm, ...
    "seconds", toc(wall));
out.perturbation_energy = perturbation_energy;
out.interface_samples = interface_samples;
out.divergence_star = div_star;
out.curve = curve;
out.record = record;
out.difference = difference;
end

function D = fluid_data(k, curve_0, amplitude, t, dt, uu, pp, physics)
[fp, ~, ~, pd] = translating_circle_problems(k, curve_0, amplitude, t, dt, uu, pp, ["0"; "0"], physics{:});
D = struct("fp", fp, "pd", pd);
end

function v = component(h, x, y, t, i)
% component i of the stacked handle h(x, y, t) at points of any shape
s = h(x(:), y(:), t);
n = numel(x);
v = reshape(s((i - 1)*n + (1:n)), size(x));
end

function [A3, A4, Ab] = wall_rows(geo, A3, A4, Ab, wall, g, d1)
% the trace rows of the wall faces (wall: 3 x n_el) replaced by the face L2
% projections of the data g (face points x n_el x [u_x, u_y, v_x, v_y]):
% <lambda - g, mu>_F = 0 for u_hat and v_hat
nq = size(geo.wa, 1)/3;
for e = 1:3
    els = find(wall(e, :));
    if isempty(els), continue, end
    nodes = (e - 1)*nq + (1:nq);
    Le = geo.Lf(nodes, (e - 1)*d1 + (1:d1), els);
    W = reshape(geo.wa(nodes, els), nq, 1, []);
    M = pagemtimes(Le, "transpose", W.*Le, "none");
    for f = 1:4
        r = (f - 1)*3*d1 + (e - 1)*d1 + (1:d1);
        A3(r, :, els) = 0;
        A4(r, :, els) = 0;
        A4(r, r, els) = M;
        Ab(r, 1, els) = -pagemtimes(Le, "transpose", W.*reshape(g(nodes, els, f), nq, 1, []), "none");
    end
end
end

function conv = convection_source(uu, rho)
% handles (x, y, t) of the two components of rho (u.grad) u, u = (-d_y,
% d_x) of the stream function uu
syms x y t real
psi = str2sym(uu);
u = [-diff(psi, y); diff(psi, x)];
c = rho*jacobian(u, [x, y])*u;
conv = arrayfun(@(e) matlabFunction(e, Vars = [x, y, t]), c, UniformOutput = false);
end

function D = with_convection(D, conv, t)
% the fluid source of D plus the convective term at time t
for c = 1:2
    f = D.fp.Data.source_term{c};
    g = conv{c};
    D.fp.Data.source_term{c} = @(x, y) f(x, y) + g(x, y, t) + 0*x;
end
end

function n = curve_normal(curve, x, y)
% the curve's unit normal (outward from the solid) at the closest points
% of (x, y), as {n_x, n_y} of the shape of x
t = curve.map_points_to_interface([reshape(x, 1, []); reshape(y, 1, [])]);
N = reshape(curve.unit_normal(reshape(t, 1, [])), 2, []);
n = {reshape(N(1, :), size(x)), reshape(N(2, :), size(x))};
end

function D = fluid_data_flat(k, y0, amplitude, t, dt, uu, pp, physics)
[fp, ~, ~, pd] = translating_flat_problems(k, y0, amplitude, t, dt, uu, pp, ["0"; "0"], physics{:});
D = struct("fp", fp, "pd", pd);
end

function v = traction(ex, sigma_F, nS, x, y, t, i)
% component i of sigma_F n_F + sigma_S n_S, n_S the solid's outward unit
% normal on Sigma, n_F = -n_S; sigma_S rows (11, 21, 12, 22)
v = 0*x;
for j = 1:2
    sS = component(ex.sigma, x, y, t, i + 2*(j - 1));
    v = v + (sS - (sigma_F{i, j}(x, y) + 0*x)).*nS{j};
end
end

function blocks = declarative_blocks(A1, A2, Af, A3, A4, Ab, p, d2, d1, h, n_el)
% the helper's element system (local [sigma; u; gamma; v], rows [C; K; S;
% M], A1 X = A2 Lam + Af; trace rows [Fu; Fv], A3 X + A4 Lam + Ab = 0,
% Lam = [u_hat (x faces, y faces); v_hat]) as declarative blocks on the
% spaces of interface_problem's Eulerian solid: local ids [sigma, w,
% gamma, u] (tests [C, M, S, K]), global ids 5 = w_hat (tests Fv), 6 =
% u_hat (tests Fu), with w = h v
perm = [1:p, p + 3*d2 + (1:2*d2), p + 2*d2 + (1:d2), p + (1:2*d2)];
iv = p + 3*d2 + (1:2*d2);
A1(:, iv, :) = A1(:, iv, :)/h;
A3(:, iv, :) = A3(:, iv, :)/h;
A2(:, 6*d1 + 1:end, :) = A2(:, 6*d1 + 1:end, :)/h;
A4(:, 6*d1 + 1:end, :) = A4(:, 6*d1 + 1:end, :)/h;
A1 = A1(perm, perm, :);
A2 = A2(perm, :, :);
Af = Af(perm, :, :);
A3 = A3(:, perm, :);
sizes = [p, 2*d2, d2, 2*d2];
first = cumsum([0, sizes(1:end - 1)]);
loc = @(i) first(i) + (1:sizes(i));
volume = @(A) reshape(A, size(A, 1), size(A, 2), 1, 1, n_el);
% columns [x faces 1..3, y faces 1..3] -> 2 d1 columns per face, face in dim 7
hybrid = @(A) permute(reshape(permute(reshape(A, size(A, 1), d1, 3, 2, n_el), ...
    [1 2 4 3 5]), size(A, 1), 2*d1, 3, n_el), [1 2 5 6 4 7 3]);
% rows [x faces 1..3, y faces 1..3] -> 2 d1 rows per face, face in dim 6
face_rows = @(A) reshape(permute(reshape(A, d1, 3, 2, size(A, 2), n_el), [1 3 4 5 2]), ...
    2*d1, size(A, 2), 1, 1, n_el, 3);
hat = {6*d1 + (1:6*d1), 1:6*d1};   % Lam columns / trace rows of ids 5 (v_hat), 6 (u_hat)
lhs = [];
flux = [];
rhs = [];
for i = 1:4
    for j = 1:4
        B = A1(loc(i), loc(j), :);
        if any(B(:)), lhs = add(lhs, IndexedArray(volume(B), [j, i])); end
    end
    for g = 1:2
        B = -A2(loc(i), hat{g}, :);
        if any(B(:)), lhs = add(lhs, IndexedArray(hybrid(B), [4 + g, i])); end
    end
    rhs = add(rhs, IndexedArray(volume(Af(loc(i), :, :)), [-1, i]));
end
for g = 1:2
    rows = hat{g};
    for j = 1:4
        B = A3(rows, loc(j), :);
        if any(B(:)), flux = add(flux, IndexedArray(face_rows(B), [j, 4 + g])); end
    end
    for gg = 1:2
        B = A4(rows, hat{gg}, :);
        if any(B(:)), flux = add(flux, IndexedArray(face_diagonal(B, d1, n_el), [4 + gg, 4 + g])); end
    end
    B = -Ab(rows, :, :);
    if any(B(:)), rhs = add(rhs, IndexedArray(face_rows(B), [-1, 4 + g])); end
end
blocks = struct("lhs", lhs, "flux", flux, "rhs", rhs);
end

function s = add(s, a)
if isempty(s), s = a; else, s = s + a; end
end

function F = face_diagonal(A, d1, n_el)
% a face-diagonal trace block (rows and columns [x faces 1..3, y faces
% 1..3]) as 2 d1 x 2 d1 per face, face in dim 6
F = zeros(2*d1, 2*d1, 1, 1, n_el, 3);
for e = 1:3
    r = [(e - 1)*d1 + (1:d1), 3*d1 + (e - 1)*d1 + (1:d1)];
    F(:, :, 1, 1, :, e) = reshape(A(r, r, :), 2*d1, 2*d1, 1, 1, n_el);
end
end

function U = project_handles(mesh, k, h)
% L2 projection onto P_k^2 of the cell of handles h{1}, h{2} of (x, y)
[rs, w] = get_2d_reference_quadrature(k + 3);
w = w(:)/sum(w);
n_el = size(mesh.elements, 2);
xy = mesh.get_2d_physical(repmat(rs, 1, 1, n_el), 1:n_el);
x = reshape(xy(1, :, :), numel(w), []);
y = reshape(xy(2, :, :), numel(w), []);
P = dubiner_basis(k, rs);
fit = (P.'*(w.*P)) \ (P.'.*w.');
U = [fit*(h{1}(x, y) + 0*x); fit*(h{2}(x, y) + 0*x)];
end

function [e, nrm] = vector_error(mesh, k, V, v)
% L2 error and exact norm of the Dubiner field V against the handles v{1}, v{2}
[rs, w] = get_2d_reference_quadrature(k + 3);
w = w(:)/sum(w);
n_el = size(mesh.elements, 2);
xy = mesh.get_2d_physical(repmat(rs, 1, 1, n_el), 1:n_el);
x = reshape(xy(1, :, :), numel(w), []);
y = reshape(xy(2, :, :), numel(w), []);
D = dubiner_basis(k, rs);
d2 = size(D, 2);
r1 = D*V(1:d2, :) - (v{1}(x, y) + 0*x);
r2 = D*V(d2 + (1:d2), :) - (v{2}(x, y) + 0*x);
e = sqrt(sum(sum(w.*(r1.^2 + r2.^2), 1).*mesh.volume(:).'));
nrm = sqrt(sum(sum(w.*((v{1}(x, y) + 0*x).^2 + (v{2}(x, y) + 0*x).^2), 1).*mesh.volume(:).'));
end

function [fluid, solid, keys] = carve_both(background, curve)
% the two sides of the background split by curve, interface vertices
% snapped onto it, and the keys of the unsnapped carve (element centroids)
[fluid, solid] = background.split(curve, max_normal_deviation = 90);
centroid = @(m) squeeze(mean(reshape(m.coordinates(:, m.elements), 2, 3, []), 2)).';
keys = struct("fluid", centroid(fluid), "solid", centroid(solid));
fluid = snap_to_curve(fluid, curve);
solid = snap_to_curve(solid, curve, normal_sign = -1);
end

function E = remap_fluid_state(D, old_mesh, new_mesh, k, old_keys, new_keys)
% velocity_stress_fsi_march's remap: an element present at both steps with
% the same shape keeps its polynomial; any other gets the L2 projection of
% its own old element's (changed shape) or the nearest old element's
% (switched on) polynomial
centroid = @(m) squeeze(mean(reshape(m.coordinates(:, m.elements), 2, 3, []), 2)).';
[present, from] = ismembertol(new_keys, old_keys, 1e-10, "ByRows", true);
c_old = centroid(old_mesh);
c_new = centroid(new_mesh);
same_shape = present;
same_shape(present) = vecnorm(c_new(present, :) - c_old(from(present), :), 2, 2) < 1e-14;
kept = find(same_shape);
added = find(~same_shape);
d2 = nchoosek(k + 2, 2);
[rs, w] = get_2d_reference_quadrature(2*k + 2);
P = dubiner_basis(k, rs);
fit = (P.'*(w(:).*P)) \ (P.'.*w(:).');
E = zeros(size(D, 1), size(new_mesh.elements, 2));
E(:, kept) = D(:, from(kept));
for j = 1:numel(added)
    if present(added(j))
        source = from(added(j));
    else
        [~, source] = min(vecnorm(c_old - c_new(added(j), :), 2, 2));
    end
    xy = new_mesh.get_2d_physical(rs, added(j));
    Q = dubiner_basis(k, old_mesh.get_2d_local(xy, source));
    for c = 0:size(D, 1)/d2 - 1
        rows = c*d2 + (1:d2);
        E(rows, added(j)) = fit*(Q*D(rows, source));
    end
end
end
