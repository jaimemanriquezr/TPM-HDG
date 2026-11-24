function out = velocity_stress_fsi_march(k, background, amplitude, dt, n_steps, uu_fluid, p_fluid, u_solid, options)
%VELOCITY_STRESS_FSI_MARCH  FSI march with a velocity-stress Eulerian solid, both sides carved.
%
%   out = velocity_stress_fsi_march(k, background, A, dt, n_steps, uu_fluid, p_fluid, u_solid, ...
%       scheme = "midpoint" | "bdf2", geometry = "flat", rotation = 0, fields = "standard", ...)
%
%   Sigma(t) moves rigidly, phi_t(X) = c(t) + Q(theta(t)) (X - c0), with
%   c(t) = c0 + A sin(2 pi t) e and theta(t) = rotation sin(2 pi t):
%   geometry = "flat", Sigma = {y = y0 + A sin(2 pi t)} (e = e_y, fluid
%   above, translating_flat_problems; rotation must be 0), or "circle" of
%   the given radius and centre c0 = center (e = e_x, fluid outside,
%   translating_circle_problems; the rotation changes the solid's
%   transport and data, not the carve). At every t_n both meshes are the
%   two sides of SimplicialGrid.split of the background (cut elements
%   discarded), their interface vertices moved onto Sigma(t_n) when snap =
%   true (snap_to_curve, normal_sign = -1 for the solid).
%   Fluid: Stokes, PEERS stress, velocity and pressure, the velocity
%   remapped across carves as in eulerian_solid_march (an element present
%   at both steps keeps its polynomial, one that switches on or changes
%   shape gets the L2 projection of the nearest old element's polynomial).
%   Solid: velocity v (relative to the rigid motion) and stress sigma,
%     rho D v - div sigma = f,  A D sigma - eps(v) = 0,
%   D the co-rotational material derivative along the rigid motion, A the
%   compliance; its state is (sigma, v) only, carried to the new carve by
%   the exact transport of velocity_stress_solid_march (the old carve moved
%   rigidly, v -> R v, sigma -> R sigma R' and its A-orthogonal weakly
%   symmetric projection, R = Q(theta_n+1 - theta_n)). extension_degree
%   = q (default [], the full polynomial) limits the old polynomial
%   extended into the part of a new element the moved old carve does not
%   cover to its degree-q Taylor polynomial (translated_field_between),
%   for v and sigma at every transported level.
%   Transmission on Sigma(t_n+1): u_F - v_S = T2, sigma_F n_F + sigma_S n_S
%   = T1 (a physical solution has T2 = the rigid velocity).
%   One coupled solve per step on the carve of t_n+1, with c = ct/h:
%     (c A sigma, tau) + (v, div tau) + (gamma, tau) - <v_hat, tau n> = (c A Z_sigma, tau)
%     (c rho v, w) - (div sigma, w) + <S (v - v_hat), w> = (c rho Z_v, w) + (f, w)
%     (sigma, eta) = 0,
%   the fluid's step and the TPM rows T1 and T2, both sides transferred
%   along their paths, the solid trace of v being v_hat + int_path grad v
%   . m with grad v = c A (sigma - Z_sigma) + gamma. It is a
%   LinearElasticityProblem in w = h v (density ct rho, Young modulus
%   E/ct, Newmark beta = 1, time step h, history displacement h Z_v,
%   Data.constitutive_source = ct (A Z_sigma, tau)) coupled by
%   FluidStructureCoupling.interface_problem under scheme "onestep" with
%   the solid state stress Z_sigma (the T2 history Ext(Z_sigma) is the
%   constitutive-source part of the path transfer of w) and T2 unscaled.
%   scheme = "midpoint": Crank-Nicolson on both sides with the unknowns at
%   the midpoint, h = dt/2, ct = 1, Z the transported state; the fluid a
%   backward-Euler step of length dt/2 from the remapped u_F^n; data at
%   the midpoint (fluid: the average of t_n and t_n+1 at fixed x; solid:
%   (R X(foot(x), t_n) + X(x, t_n+1))/2 along the rigid motion, R . R'
%   for the stress; T1, T2 from those averages with the normals of
%   Sigma(t_n+1)); then X^n+1 = 2 X' - X~^n for u_F, sigma_S and v_S.
%   scheme = "bdf2": the first step is a midpoint step, then BDF2 on both
%   sides: h = dt, ct = 3/2, Z = (4 X~_n - X~_n-1)/3 with X~_n-1 the
%   previous transported level transported once more; the fluid the theta
%   = 2/3 step from (4 u^n - u^n-1)/3 with zero acceleration; data at
%   t_n+1.
%   fields = "standard": body-frame displacement U(X, t) = u_S(phi_t(X),
%   t); v = Q d_t U, sigma = Q C eps_X(U) Q' at X = phi_t^-1(x) (with
%   rotation = 0, v = D_t u_S and sigma = C eps(u_S), the helpers' solid
%   data). fields = "patch": V(X) of degree k and Sigma(X, t) = S0(X) +
%   t C eps_X(V) in the body frame. The flat solid's walls are Dirichlet
%   with h v_exact.
%   solid_tau: the solid stabilization in velocity form, tau_S (v - v_hat)
%   on the element boundaries; the solve's unknown is w = h v, so the
%   solid mesh carries tau = solid_tau/h in each solve (the T1 row's solid
%   pair included). The fluid's tau is the mesh default.
%   initial_perturbation = [seed, amplitude]: adds to the initial (sigma,
%   v) a random weakly symmetric pair of energy amplitude^2 times the
%   initial solid energy.
%   start_step, start_state = out.state of a run to step start_step:
%   continues that run; the records then cover the steps after it.
%   record_state = true: out.record(i) = struct(step, sigma, v,
%   fluid_velocity), the solid and fluid dofs after each step of the run.
%   reference = such a record of another run on the same carves (it may
%   cover more steps; entries are matched by step): out.difference (per
%   step) = [1/2 (rho |dv|^2 + (A dsigma, dsigma)), rho_F/2 |du_F|^2] of
%   this run's state minus the reference's.
%
%   out: eF = [sigma_F, u_F] and e_sigma_S, e_v_S, L2 errors at t =
%   n_steps dt, n_sigma_F, n_u_F, n_sigma_S, n_v_S the exact norms there;
%   step_errors, per step [sigma_F, u_F, sigma_S, v_S, sigma_F_extrapolated]
%   (under a midpoint step sigma_F is the midpoint stress against the
%   averaged exact stress, and sigma_F_extrapolated 2 sigma_F' - sigma_F^n
%   (remapped) against sigma_F(t_n+1), a diagnostic not fed back; under
%   a BDF2 step both are the solved stress at t_n+1); energy (the solid
%   energy 1/2 (rho |v|^2 + (A sigma, sigma)) at the start and after each
%   step) and energy_transported (after each transport); uncovered (per
%   step the largest fraction of a new solid element and the total area
%   not covered by the moved old carve); history, per step [t, max|u_F|,
%   max|v_S|, fluid on, fluid off]; paths, the transfer path lengths [min,
%   max]; state (the restart state); record and difference (see above);
%   fluid_mesh, solid_mesh; h, the largest face; seconds and seconds_data
%   (symbolic data).
%
%See also: velocity_stress_solid_march, eulerian_solid_march,
%          FluidStructureCoupling.interface_problem
arguments
    k (1,1) double
    background SimplicialGrid
    amplitude (1,1) double
    dt (1,1) double
    n_steps (1,1) double
    uu_fluid (1,1) string
    p_fluid (1,1) string
    u_solid (2,1) string
    options.scheme (1,1) string {mustBeMember(options.scheme, ["midpoint", "bdf2"])} = "bdf2"
    options.geometry (1,1) string {mustBeMember(options.geometry, ["flat", "circle"])} = "flat"
    options.y0 (1,1) double = 0.45
    options.center (2,1) double = [0.5; 0.5]
    options.radius (1,1) double = 0.2
    options.rotation (1,1) double = 0
    options.snap (1,1) logical = true
    options.rho_f = 1; options.mu = 1; options.rho_s = 1; options.E = 10; options.nu = 0.3;
    options.fields (1,1) string {mustBeMember(options.fields, ["standard", "patch"])} = "standard"
    options.swap_transmission (1,1) logical = true
    options.n_quadpoints double {mustBeScalarOrEmpty} = []
    % per step errors printed
    options.verbose (1,1) logical = false
    options.initial_perturbation double = []
    options.start_step (1,1) double {mustBeInteger, mustBeNonnegative} = 0
    options.start_state struct = struct.empty
    options.solid_tau (1,1) double {mustBePositive} = 1
    options.record_state (1,1) logical = false
    options.reference struct = struct.empty
    options.extension_degree double {mustBeScalarOrEmpty} = []
end
wall = tic;
flat = options.geometry == "flat";
assert(~flat || options.rotation == 0, "velocity_stress_fsi_march:rotation", ...
    "The flat interface only translates (rotation must be 0).");
physics = {options.rho_f, options.mu, options.rho_s, options.E, options.nu};
if flat
    data = @(t) translating_flat_problems(k, options.y0, amplitude, t, dt, uu_fluid, p_fluid, ...
        u_solid, physics{:});
    curve_at = @(t) FlatCurveHorizontal(options.y0 + amplitude*sin(2*pi*t));
    direction = [0; 1];
    c0 = [0.5; options.y0];
    normal_S = @(x, y, t) {0*x; 1 + 0*x};
else
    curve_0 = CircularCurve(options.center, options.radius);
    data = @(t) translating_circle_problems(k, curve_0, amplitude, t, dt, uu_fluid, p_fluid, ...
        u_solid, physics{:});
    curve_at = @(t) CircularCurve(curve_0.center + [amplitude*sin(2*pi*t); 0], curve_0.radius);
    direction = [1; 0];
    c0 = options.center;
    normal_S = @(x, y, t) {(x - c0(1) - amplitude*sin(2*pi*t))/options.radius; ...
        (y - c0(2))/options.radius};
end
cen = @(t) c0 + amplitude*sin(2*pi*t)*direction;
ang = @(t) options.rotation*sin(2*pi*t);
rot = @(a) [cos(a), -sin(a); sin(a), cos(a)];
moving = amplitude ~= 0 || options.rotation ~= 0;
rho = options.rho_s;
[lambda, mu_s] = LinearElasticityProblem.get_lame_parameters(options.E, options.nu);
c1 = 1/(2*mu_s);
c2 = lambda/(4*mu_s*(lambda + mu_s));
ex = exact_fields(k, options, u_solid, c0, direction, amplitude, lambda, mu_s);
sp0 = LinearElasticityProblem(k);   % only the object translated_field_between needs
carve = @(curve) carve_both(background, curve, options.snap);
quad = {"n_quadpoints", options.n_quadpoints, "n_quadpoints_path", options.n_quadpoints + 1};
if isempty(options.n_quadpoints), quad = {"n_quadpoints", k + 2, "n_quadpoints_path", k + 3}; end
seconds_data = 0;
    function D = load_data(t)
        % the helper's fluid problem, interface sources and exact fluid fields at t
        t_data = tic;
        [fp_t, sp_t, src_t, pd_t] = data(t);
        D = struct("fp", fp_t, "sp", sp_t, "src", {src_t}, "pd", pd_t);
        seconds_data = seconds_data + toc(t_data);
    end

n0 = options.start_step;
[fm, sm, keys] = carve(curve_at(n0*dt));
ops = stress_operators(sm, k, c1, c2);
if n0 == 0
    D0 = load_data(0);
    fp = D0.fp;
    fp.Data.time_scheme = "euler";
    sp_h = D0.sp;
    sp_h.Data.time_scheme = "newmark";
    initial = FluidStructureCoupling.initial_state(fp, sp_h, fm, sm, curve_at(0), D0.src{:}, ...
        "Order", k, quad{:});
    UF = initial.fluid.velocity;
    SF = project_peers(fm, k, @(x, y) evaluate(D0.pd.sigma_F, x, y));
    V = project_vector(sm, k, @(x, y) ex.v(x, y, 0));
    S = project_stress(sm, k, ops, c1, c2, @(x, y) ex.sigma(x, y, 0));
    if ~isempty(options.initial_perturbation)
        stream = RandStream("mt19937ar", Seed = options.initial_perturbation(1));
        [dS, dV] = random_state(stream, sm, k, ops, rho);
        a = options.initial_perturbation(2)*sqrt(energy(sm, k, ops, rho, S, V));
        S = S + a*dS;
        V = V + a*dV;
    end
    [v_previous, previous] = deal([]);
else
    s = options.start_state;
    assert(~isempty(s) && s.step == n0, "velocity_stress_fsi_march:start", ...
        "A restart needs the state of a run to step start_step.");
    [UF, SF, v_previous, S, V, previous] = deal(s.fluid_velocity, s.fluid_stress, ...
        s.fluid_velocity_previous, s.sigma, s.v, s.previous);
    D0 = struct.empty;
    if options.scheme == "midpoint", D0 = load_data(n0*dt); end
end

m = n_steps - n0;
E = nan(m + 1, 1);
E_t = nan(m, 1);
E(1) = energy(sm, k, ops, rho, S, V);
uncovered = zeros(m, 2);
step_errors = nan(m, 5);
history = zeros(m, 5);
paths = [Inf, 0];
record = struct("step", {}, "sigma", {}, "v", {}, "fluid_velocity", {});
difference = nan(m, 2);
for n = n0 + 1:n_steps
    t0 = (n - 1)*dt;
    t1 = n*dt;
    midpoint = options.scheme == "midpoint" || n == 1;
    D1 = load_data(t1);
    curve = curve_at(t1);
    [fm_old, sm_old, keys_old] = deal(fm, sm, keys);
    [fm, sm, keys] = carve(curve);
    % fluid: velocity (and the previous velocity and the diagnostic stress)
    % remapped to the new carve
    [fluid, switched] = remap_fluid_state(struct("velocity", UF), fm_old, fm, k, keys_old.fluid, ...
        keys.fluid, "velocity");
    UF = fluid.velocity;
    SF = remap_fluid_stress(SF, fm_old, fm, k, keys_old.fluid, keys.fluid);
    if ~isempty(v_previous)
        v_previous = remap_fluid_state(struct("velocity", v_previous), fm_old, fm, k, ...
            keys_old.fluid, keys.fluid, "velocity").velocity;
    end
    % solid: (sigma, v) transported to the new carve
    R = rot(ang(t1) - ang(t0));
    [c_old, c_new] = deal(cen(t0), cen(t1));
    if moving
        ops = stress_operators(sm, k, c1, c2);
        [Vt, St, uncovered(n - n0, :)] = transport(sp0, sm_old, sm, k, V, S, c_old, c_new, R, ops, ...
            c1, c2, options.extension_degree);
        if ~isempty(previous)
            [previous.V, previous.S] = transport(sp0, sm_old, sm, k, previous.V, previous.S, ...
                c_old, c_new, R, ops, c1, c2, options.extension_degree);
        end
    else
        [Vt, St] = deal(V, S);
    end
    E_t(n - n0) = energy(sm, k, ops, rho, St, Vt);

    fp = D1.fp;
    pd = D1.pd;
    if midpoint
        h = dt/2;
        ct = 1;
        fp.Data.time_scheme = "euler";
        fp.Data.theta = 1;
        fp.Data.source_term = average(D0.fp.Data.source_term, D1.fp.Data.source_term);
        fp.Data.datum_dirichlet = average(D0.fp.Data.datum_dirichlet, D1.fp.Data.datum_dirichlet);
        fluid_history = UF;
        sigma_F = average(D0.pd.sigma_F, pd.sigma_F);
        u_F = average(D0.pd.u_F, pd.u_F);
        [Zv, Zs] = deal(Vt, St);
        f = averaged(ex.f, t0, t1, R, c_old, c_new);
        g = averaged(ex.v, t0, t1, R, c_old, c_new);
        sigma_S = averaged_tensor(ex.sigma, t0, t1, R, c_old, c_new);
    else
        h = dt;
        ct = 3/2;
        fp.Data.time_scheme = "cranknicolson";
        fp.Data.theta = 2/3;
        fluid_history = (4*UF - v_previous)/3;
        [sigma_F, u_F] = deal(pd.sigma_F, pd.u_F);
        [Zv, Zs] = deal((4*Vt - previous.V)/3, (4*St - previous.S)/3);
        f = at_time(ex.f, t1);
        g = at_time(ex.v, t1);
        sigma_S = at_time_tensor(ex.sigma, t1);
    end
    fp.Data.time_step = h;
    sp = LinearElasticityProblem(k, density = ct*rho, time_step = h, time_scheme = "newmark", ...
        newmark_beta = 1, newmark_gamma = 0.5, young_modulus = options.E/ct, ...
        poisson_ratio = options.nu, source_term = f, ...
        datum_dirichlet = {@(x, y) h*g{1}(x, y); @(x, y) h*g{2}(x, y)});
    sp.Data.constitutive_source = ct*reshape(pagemtimes(ops.T, reshape(Zs, [], 1, size(Zs, 2))), ...
        size(Zs));
    solve_state = struct("fluid", struct("velocity", fluid_history, ...
        "acceleration", zeros(size(fluid_history))), "solid", struct("displacement", h*Zv, ...
        "velocity", 0, "acceleration", 0, "stress", Zs, "rotation", 0));
    % T1 = sigma_F n_F + sigma_S n_S (n_F = -n_S), T2 = u_F - v_S
    T1 = cell(2, 1);
    T2 = cell(2, 1);
    for i = 1:2
        T1{i} = @(x, y) traction(sigma_S, sigma_F, normal_S, x, y, t1, i);
        T2{i} = @(x, y) u_F{i}(x, y) - g{i}(x, y);
    end
    sm_solve = sm.add_stabilization(options.solid_tau/h);
    coupled = FluidStructureCoupling.interface_problem(fp, sp, fm, sm_solve, curve, T1, T2, "onestep", ...
        "Order", k, "state", solve_state, "patch_weight", "jacobian", quad{:}, ...
        "swap_transmission", options.swap_transmission);
    [fs, ss, gs] = solve(coupled, k);
    [Yf, Ysf] = deal(fs.Unknowns.field(2).dof, fs.Unknowns.field(1).dof);
    [Ys, Yv] = deal(ss.Unknowns.field(1).dof, ss.Unknowns.field(2).dof/h);
    v_previous = UF;
    if midpoint
        [UF, SF, S, V] = deal(2*Yf - UF, 2*Ysf - SF, 2*Ys - St, 2*Yv - Vt);
    else
        [UF, SF, S, V] = deal(Yf, Ysf, Ys, Yv);
    end
    if options.scheme == "bdf2"
        previous = struct("V", Vt, "S", St);
    end
    D0 = D1;

    E(n - n0 + 1) = energy(sm, k, ops, rho, S, V);
    exact_sigma_F = @(x, y) evaluate(pd.sigma_F, x, y);
    [eSF, nSF] = tensor_error(fm, k, Ysf, @(x, y) evaluate(sigma_F, x, y));
    [eSFx, nSFx] = tensor_error(fm, k, SF, exact_sigma_F);
    [eUF, nUF] = vector_error(fm, k, UF, @(x, y) evaluate(pd.u_F, x, y));
    [eSS, eVS, nSS, nVS] = errors(sm, k, S, V, ex, t1);
    step_errors(n - n0, :) = [eSF, eUF, eSS, eVS, eSFx];
    history(n - n0, :) = [t1, max(abs(UF), [], "all"), max(abs(V), [], "all"), switched];
    if options.record_state
        record(n - n0) = struct("step", n, "sigma", S, "v", V, "fluid_velocity", UF);
    end
    if ~isempty(options.reference)
        r = options.reference([options.reference.step] == n);
        difference(n - n0, :) = [energy(sm, k, ops, rho, S - r.sigma, V - r.v), ...
            fluid_energy(fm, k, options.rho_f, UF - r.fluid_velocity)];
    end
    for part = {gs.partition_1, gs.partition_2}
        l = vecnorm(part{1}.quadrature.interface - part{1}.quadrature.physical, 2, 1);
        paths = [min(paths(1), min(l(:))), max(paths(2), max(l(:)))];
    end
    if options.verbose
        fprintf("  step %d t %.4f %s  sF %.3e  uF %.3e  sS %.3e  vS %.3e  sFx %.3e  on/off F %d/%d\n", ...
            n, t1, ternary(midpoint, "mid ", "bdf2"), step_errors(n - n0, :), switched);
    end
end
state = struct("step", n_steps, "fluid_velocity", UF, "fluid_velocity_previous", v_previous, ...
    "fluid_stress", SF, "sigma", S, "v", V, "previous", previous);
if m == 0
    [eSF, eUF, eSS, eVS, nSF, nUF, nSS, nVS] = deal(NaN);
end
out = struct("eF", [eSF, eUF], "e_sigma_S", eSS, "e_v_S", eVS, "n_sigma_F", nSF, "n_u_F", nUF, ...
    "n_sigma_S", nSS, "n_v_S", nVS, "n_sigma_F_extrapolated", nSFx, "step_errors", step_errors, ...
    "energy", E, "energy_transported", E_t, "uncovered", uncovered, "history", history, ...
    "paths", paths, "state", state, "record", record, "difference", difference, "fluid_mesh", fm, "solid_mesh", sm, ...
    "h", max([fm.area(:); sm.area(:)]), "seconds", toc(wall), "seconds_data", seconds_data);
end

function r = ternary(c, a, b)
if c, r = a; else, r = b; end
end

function v = traction(sigma_S, sigma_F, normal_S, x, y, t, i)
% component i of sigma_F n_F + sigma_S n_S, n_F = -n_S
n = normal_S(x, y, t);
v = 0*x;
for j = 1:2
    v = v + (sigma_S{i, j}(x, y) - sigma_F{i, j}(x, y)).*n{j};
end
end

function c = evaluate(handles, x, y)
c = cellfun(@(f) f(x, y) + 0*x, handles, UniformOutput = false);
end

function a = average(h0, h1)
% the cell of handles (h0 + h1)/2
a = cell(size(h0));
for i = 1:numel(h0)
    [f, g] = deal(h0{i}, h1{i});
    a{i} = @(x, y) (f(x, y) + g(x, y))/2;
end
end

function [fluid, solid, keys] = carve_both(background, curve, snap)
% the two sides of the background split by curve and the keys of the
% unsnapped carve (element centroids), as eulerian_solid_march's carve_both
split = @(grid) grid.split(curve);
if snap
    split = @(grid) grid.split(curve, max_normal_deviation = 90);
end
[fluid, solid] = split(background);
centroid = @(m) squeeze(mean(reshape(m.coordinates(:, m.elements), 2, 3, []), 2)).';
keys = struct("fluid", centroid(fluid), "solid", centroid(solid));
if snap
    fluid = snap_to_curve(fluid, curve);
    solid = snap_to_curve(solid, curve, normal_sign = -1);
end
end

function [kept, from, added, source, switched] = fluid_sources(old_mesh, new_mesh, old_keys, new_keys)
% eulerian_solid_march's remap rule: an element present at both steps with
% the same shape keeps its polynomial; any other takes its own old element
% (changed shape) or the old element with the nearest centroid (switched on)
centroid = @(m) squeeze(mean(reshape(m.coordinates(:, m.elements), 2, 3, []), 2)).';
[present, from] = ismembertol(new_keys, old_keys, 1e-10, "ByRows", true);
switched = [nnz(~present), size(old_mesh.elements, 2) - nnz(present)];
c_old = centroid(old_mesh);
c_new = centroid(new_mesh);
same_shape = present;
same_shape(present) = vecnorm(c_new(present, :) - c_old(from(present), :), 2, 2) < 1e-14;
kept = find(same_shape);
added = find(~same_shape);
source = zeros(size(added));
for j = 1:numel(added)
    if present(added(j))
        source(j) = from(added(j));
    else
        [~, source(j)] = min(vecnorm(c_old - c_new(added(j), :), 2, 2));
    end
end
end

function [s, switched] = remap_fluid_state(s, old_mesh, new_mesh, k, old_keys, new_keys, names)
% the fluid Dubiner volume fields names from old_mesh to new_mesh, as
% eulerian_solid_march's remap_fluid_state
[kept, from, added, source, switched] = fluid_sources(old_mesh, new_mesh, old_keys, new_keys);
d2 = nchoosek(k + 2, 2);
[rs, w] = get_2d_reference_quadrature(2*k + 2);
P = dubiner_basis(k, rs);
fit = (P.'*(w(:).*P)) \ (P.'.*w(:).');
for name = names
    D = s.(name);
    E = zeros(size(D, 1), size(new_mesh.elements, 2));
    E(:, kept) = D(:, from(kept));
    for j = 1:numel(added)
        xy = new_mesh.get_2d_physical(rs, added(j));
        rs_old = old_mesh.get_2d_local(xy, source(j));
        Q = dubiner_basis(k, rs_old);
        for c = 0:size(D, 1)/d2 - 1
            rows = c*d2 + (1:d2);
            E(rows, added(j)) = fit*(Q*D(rows, source(j)));
        end
    end
    s.(name) = E;
end
end

function S = remap_fluid_stress(S, old_mesh, new_mesh, k, old_keys, new_keys)
% the PEERS fluid stress by the same rule: kept elements keep their dofs,
% the others get the L2 projection onto the PEERS space of the source
% element's field (a polynomial of degree k + 1) evaluated on them
[kept, from, added, source] = fluid_sources(old_mesh, new_mesh, old_keys, new_keys);
E = zeros(size(S, 1), size(new_mesh.elements, 2));
E(:, kept) = S(:, from(kept));
if ~isempty(added)
    [rs, w] = rule(k + 3);
    xy = new_mesh.get_2d_physical(repmat(rs, 1, 1, numel(added)), added(:).');
    G_old = get_affine_map(old_mesh, source(:).', inverse = true);
    vals = zeros(4*numel(w), numel(added));
    for j = 1:numel(added)
        r_old = old_mesh.get_2d_local(xy(:, :, j), source(j));
        P = PEERS_basis(k, r_old(1, :), r_old(2, :), G_old(:, :, j));
        Pp = reshape(permute(P, [1 3 4 2]), 4*numel(w), []);
        vals(:, j) = Pp*S(:, source(j));
    end
    E(:, added) = peers_fit(new_mesh, k, vals, rs, w, added);
end
S = E;
end

function S = peers_fit(mesh, k, vals, rs, w, elements)
% L2 projection onto the PEERS space of the given elements of the values
% ((nq*4) x n) at the points rs
G = get_affine_map(mesh, elements(:).', inverse = true);
P = PEERS_basis(k, rs(1, :), rs(2, :), G);
nq = size(P, 1);
Pp = reshape(permute(P, [1 3 4 2 5]), nq*4, size(P, 2), []);
W = repmat(w, 4, 1);
M = pagemtimes(Pp, "transpose", W.*Pp, "none");
b = pagemtimes(Pp, "transpose", reshape(W.*vals, [], 1, numel(elements)), "none");
S = reshape(pagemldivide(M, b), size(P, 2), []);
end

function S = project_peers(mesh, k, sigma)
% L2 projection onto the PEERS space of the handle sigma(x, y) {2 x 2}
[rs, w] = rule(k + 3);
n_el = size(mesh.elements, 2);
xy = mesh.get_2d_physical(repmat(rs, 1, 1, n_el), 1:n_el);
s = sigma(reshape(xy(1, :, :), numel(w), []), reshape(xy(2, :, :), numel(w), []));
S = peers_fit(mesh, k, [s{1, 1}; s{2, 1}; s{1, 2}; s{2, 2}], rs, w, 1:n_el);
end

function [e, nrm] = tensor_error(mesh, k, S, sigma)
% L2 error and exact norm of the PEERS field S against sigma(x, y) {2 x 2}
[rs, w] = rule(k + 3);
n_el = size(mesh.elements, 2);
xy = mesh.get_2d_physical(repmat(rs, 1, 1, n_el), 1:n_el);
s = sigma(reshape(xy(1, :, :), numel(w), []), reshape(xy(2, :, :), numel(w), []));
s = [s{1, 1}; s{2, 1}; s{1, 2}; s{2, 2}];
sh = stress_values(mesh, k, S, rs);
W4 = repmat(w, 4, 1);
e = sqrt(sum(sum(W4.*(sh - s).^2, 1).*mesh.volume));
nrm = sqrt(sum(sum(W4.*s.^2, 1).*mesh.volume));
end

function [e, nrm] = vector_error(mesh, k, V, v)
% L2 error and exact norm of the Dubiner field V against v(x, y) {2}
[rs, w] = rule(k + 3);
n_el = size(mesh.elements, 2);
xy = mesh.get_2d_physical(repmat(rs, 1, 1, n_el), 1:n_el);
v = v(reshape(xy(1, :, :), numel(w), []), reshape(xy(2, :, :), numel(w), []));
D = dubiner_basis(k, rs);
d2 = size(D, 2);
vh = {D*V(1:d2, :), D*V(d2 + (1:d2), :)};
e = sqrt(sum(sum(w.*((vh{1} - v{1}).^2 + (vh{2} - v{2}).^2), 1).*mesh.volume));
nrm = sqrt(sum(sum(w.*(v{1}.^2 + v{2}.^2), 1).*mesh.volume));
end

function ex = exact_fields(k, options, u_solid, c0, direction, amplitude, lambda, mu)
% handles v(x, y, t) {2}, sigma(x, y, t) {2 x 2}, f(x, y, t) {2} of the
% solid (velocity_stress_solid_march's construction)
syms x y t X Y real
c = c0 + amplitude*sin(2*pi*t)*direction;
th = options.rotation*sin(2*pi*t);
Q = [cos(th), -sin(th); sin(th), cos(th)];
C_eps = @(U) 2*mu*sym_eps(U, X, Y) + lambda*trace(sym_eps(U, X, Y))*eye(2);
switch options.fields
    case "standard"
        u = str2sym(u_solid(:));
        phi = c + Q*([X; Y] - c0);
        U = subs(u, [x; y], phi);
        Vb = diff(U, t);
        Sb = C_eps(U);
    case "patch"
        Vb = sym(zeros(2, 1));
        S0 = sym(zeros(2, 2));
        i = 0;
        for a = 0:k
            for b = 0:k - a
                i = i + 1;
                m = X^a*Y^b;
                Vb = Vb + m*[cos(1 + i); sin(2 + 3*i)];
                s = [cos(3 + i), sin(5 + 2*i); 0, cos(7 + 5*i)];
                S0 = S0 + m*(s + triu(s, 1).');
            end
        end
        Sb = S0 + t*C_eps(Vb);
end
Fb = options.rho_s*diff(Vb, t) - [diff(Sb(1, 1), X) + diff(Sb(1, 2), Y); ...
    diff(Sb(2, 1), X) + diff(Sb(2, 2), Y)];
Xx = c0 + Q.'*([x; y] - c);
to_x = @(e) subs(e, [X; Y], Xx);
v = to_x(Q*Vb);
sigma = to_x(Q*Sb*Q.');
f = to_x(Q*Fb);
handle = @(e) vectorized(matlabFunction(e, Vars = [x, y, t]));
hv = arrayfun(handle, v, UniformOutput = false);
hs = arrayfun(handle, sigma, UniformOutput = false);
hf = arrayfun(handle, f, UniformOutput = false);
ex.v = @(xx, yy, tt) cellfun(@(h) h(xx, yy, tt), hv, UniformOutput = false);
ex.sigma = @(xx, yy, tt) cellfun(@(h) h(xx, yy, tt), hs, UniformOutput = false);
ex.f = @(xx, yy, tt) cellfun(@(h) h(xx, yy, tt), hf, UniformOutput = false);
end

function E = sym_eps(U, X, Y)
G = [diff(U(1), X), diff(U(1), Y); diff(U(2), X), diff(U(2), Y)];
E = (G + G.')/2;
end

function h = vectorized(f)
h = @(x, y, t) f(x, y, t) + 0*x;
end

function h = at_time(F, t)
% {2} handles (x, y) of a vector field at time t
h = {@(x, y) pick(F(x, y, t), 1); @(x, y) pick(F(x, y, t), 2)};
end

function h = at_time_tensor(F, t)
% {2 x 2} handles (x, y) of a tensor field at time t
h = cell(2, 2);
for i = 1:4
    h{i} = @(x, y) pick(F(x, y, t), i);
end
end

function h = averaged(F, t0, t1, R, c_old, c_new)
% {2} handles (x, y) of (R F(foot(x), t0) + F(x, t1))/2
h = cell(2, 1);
for i = 1:2
    h{i} = @(x, y) (rotated_at_foot(F, x, y, t0, R, c_old, c_new, i) + pick(F(x, y, t1), i))/2;
end
end

function h = averaged_tensor(F, t0, t1, R, c_old, c_new)
% {2 x 2} handles (x, y) of (R F(foot(x), t0) R' + F(x, t1))/2
h = cell(2, 2);
for i = 1:2
    for j = 1:2
        h{i, j} = @(x, y) (rotated_tensor_at_foot(F, x, y, t0, R, c_old, c_new, i, j) ...
            + pick(F(x, y, t1), sub2ind([2 2], i, j)))/2;
    end
end
end

function [fx, fy] = foot(x, y, R, c_old, c_new)
fx = c_old(1) + R(1, 1)*(x - c_new(1)) + R(2, 1)*(y - c_new(2));
fy = c_old(2) + R(1, 2)*(x - c_new(1)) + R(2, 2)*(y - c_new(2));
end

function w = rotated_at_foot(F, x, y, t, R, c_old, c_new, i)
[fx, fy] = foot(x, y, R, c_old, c_new);
Ff = F(fx, fy, t);
w = R(i, 1)*Ff{1} + R(i, 2)*Ff{2};
end

function w = rotated_tensor_at_foot(F, x, y, t, R, c_old, c_new, i, j)
[fx, fy] = foot(x, y, R, c_old, c_new);
Ff = F(fx, fy, t);
w = 0*x;
for p = 1:2
    for q = 1:2
        w = w + R(i, p)*Ff{p, q}*R(j, q);
    end
end
end

function v = pick(c, i)
v = c{i};
end

function [rs, w] = rule(q)
[rs, w] = get_2d_reference_quadrature(q);
w = w(:)/sum(w);
end

function ops = stress_operators(mesh, k, c1, c2)
% per element: T = (A sigma, tau) matrix, the weak symmetry block MA and
% the saddle matrix of the A-orthogonal weakly symmetric projection
[M_P, M_tr, M_A] = PEERS_matrix(mesh, k);
ops.T = c1*M_P - c2*M_tr;
nP = size(M_P, 1);
d2 = size(M_A, 2);
ops.K = [ops.T, M_A; permute(M_A, [2 1 3]), zeros(d2, d2, size(M_A, 3))];
ops.nP = nP;
end

function S = sym_projection(ops, rhs)
% sigma in S_h with (A sigma, tau) + (xi, tau) = rhs(tau), (sigma, eta) = 0
n_el = size(rhs, 2);
d2 = size(ops.K, 1) - ops.nP;
X = pagemldivide(ops.K, cat(1, reshape(rhs, ops.nP, 1, n_el), zeros(d2, 1, n_el)));
S = reshape(X(1:ops.nP, 1, :), ops.nP, n_el);
end

function Pp = peers_points(mesh, k, rs)
% PEERS basis at the reference points rs of every element, as
% (nq*4) x nP x n_el, rows q fastest then the components (1,1), (2,1),
% (1,2), (2,2)
G = get_affine_map(mesh, inverse = true);
P = PEERS_basis(k, rs(1, :), rs(2, :), G);
[nq, nP] = size(P, 1:2);
Pp = reshape(permute(P, [1 3 4 2 5]), nq*4, nP, []);
end

function vals = stress_values(mesh, k, S, rs)
% (nq*4) x n_el values of the PEERS fields S at rs
Pp = peers_points(mesh, k, rs);
vals = reshape(pagemtimes(Pp, reshape(S, size(S, 1), 1, [])), size(Pp, 1), []);
end

function rhs = stress_moments(mesh, k, Avals, rs, w)
% (A sigma, tau_j)_K from the values of A sigma ((nq*4) x n_el) at rs
Pp = peers_points(mesh, k, rs);
W = repmat(w, 4, 1);
rhs = reshape(pagemtimes(Pp, "transpose", reshape(W.*Avals, [], 1, size(Avals, 2)), "none"), ...
    size(Pp, 2), []).*mesh.volume;
end

function Avals = compliance(vals, c1, c2)
% A sigma at the points, values (nq*4) x n_el in the component order
nq = size(vals, 1)/4;
tr = vals(1:nq, :) + vals(3*nq + (1:nq), :);
Avals = c1*vals;
Avals(1:nq, :) = Avals(1:nq, :) - c2*tr;
Avals(3*nq + (1:nq), :) = Avals(3*nq + (1:nq), :) - c2*tr;
end

function vals = rotate_tensor(vals, R)
% R sigma R' at every point
nq = size(vals, 1)/4;
s = @(i) vals((i - 1)*nq + (1:nq), :);
S = {s(1), s(3); s(2), s(4)};   % S{a, b} = sigma_ab
out = cell(2, 2);
for a = 1:2
    for b = 1:2
        out{a, b} = 0;
        for p = 1:2
            for q = 1:2
                out{a, b} = out{a, b} + R(a, p)*S{p, q}*R(b, q);
            end
        end
    end
end
vals = [out{1, 1}; out{2, 1}; out{1, 2}; out{2, 2}];
end

function S = project_stress(mesh, k, ops, c1, c2, sigma)
% A-orthogonal weakly symmetric projection of the handle sigma(x, y) {2 x 2}
[rs, w] = rule(k + 3);
xy = mesh.get_2d_physical(repmat(rs, 1, 1, size(mesh.elements, 2)), 1:size(mesh.elements, 2));
s = sigma(reshape(xy(1, :, :), numel(w), []), reshape(xy(2, :, :), numel(w), []));
vals = [s{1, 1}; s{2, 1}; s{1, 2}; s{2, 2}];
S = sym_projection(ops, stress_moments(mesh, k, compliance(vals, c1, c2), rs, w));
end

function V = project_vector(mesh, k, v)
% L2 projection onto P_k of the handle v(x, y) {2}
[rs, w] = rule(k + 3);
xy = mesh.get_2d_physical(repmat(rs, 1, 1, size(mesh.elements, 2)), 1:size(mesh.elements, 2));
s = v(reshape(xy(1, :, :), numel(w), []), reshape(xy(2, :, :), numel(w), []));
D = dubiner_basis(k, rs);
fit = (D.'*(w.*D)) \ (D.'.*w.');
V = [fit*s{1}; fit*s{2}];
end

function [S, V] = random_state(stream, mesh, k, ops, rho)
% random (sigma, v) of unit energy, sigma weakly symmetric
n_el = size(mesh.elements, 2);
d2 = nchoosek(k + 2, 2);
V = randn(stream, 2*d2, n_el);
S0 = randn(stream, ops.nP, n_el);
S = sym_projection(ops, reshape(pagemtimes(ops.T, reshape(S0, ops.nP, 1, n_el)), ops.nP, n_el));
a = 1/sqrt(energy(mesh, k, ops, rho, S, V));
S = a*S;
V = a*V;
end

function E = energy(mesh, k, ops, rho, S, V)
% 1/2 (rho |v|^2 + (A sigma, sigma))
[rs, w] = rule(k + 1);
D = dubiner_basis(k, rs);
M = D.'*(w.*D);
d2 = size(M, 1);
ev = 0;
for c = 0:1
    Vc = V(c*d2 + (1:d2), :);
    ev = ev + sum(sum(Vc.*(M*Vc), 1).*mesh.volume);
end
n_el = size(S, 2);
es = sum(pagemtimes(reshape(S, [], 1, n_el), "transpose", ...
    pagemtimes(ops.T, reshape(S, [], 1, n_el)), "none"), "all");
E = (rho*ev + es)/2;
end

function E = fluid_energy(mesh, k, rho, U)
% rho/2 |u|^2 of the Dubiner field U
[rs, w] = rule(k + 1);
D = dubiner_basis(k, rs);
M = D.'*(w.*D);
d2 = size(M, 1);
E = 0;
for c = 0:1
    Uc = U(c*d2 + (1:d2), :);
    E = E + sum(sum(Uc.*(M*Uc), 1).*mesh.volume);
end
E = rho*E/2;
end

function [Vt, St, uncovered] = transport(sp0, old, new, k, V, S, c_old, c_new, R, ops, c1, c2, q)
% (sigma, v) of old carried by the rigid step onto new
% (velocity_stress_solid_march's transport; q the extension degree)
moved = old;
X = old.coordinates;
moved.coordinates = c_new + R*(X - c_old);
nn = reshape(old.normals, 2, []);
moved.normals = reshape(R*nn, size(old.normals));
[Vt, info] = sp0.translated_field_between(moved, new, k, k, V, [0; 0], extension_degree = q);
d2 = nchoosek(k + 2, 2);
v1 = Vt(1:d2, :);
v2 = Vt(d2 + (1:d2), :);
Vt = [R(1, 1)*v1 + R(1, 2)*v2; R(2, 1)*v1 + R(2, 2)*v2];
area = new.volume(:);
uncovered = [max(info.uncovered), sum(info.uncovered(:).*area)];
[rs, w] = rule(k + 2);
Dk1 = dubiner_basis(k + 1, rs);
fit = (Dk1.'*(w.*Dk1)) \ (Dk1.'.*w.');
vals = stress_values(old, k, S, rs);
nq = numel(w);
F = zeros(4*size(fit, 1), size(S, 2));
d3 = size(fit, 1);
for i = 1:4
    F((i - 1)*d3 + (1:d3), :) = fit*vals((i - 1)*nq + (1:nq), :);
end
assert(norm(kron(eye(4), Dk1)*F - vals, "fro") <= 1e-9*max(1, norm(vals, "fro")), ...
    "velocity_stress_fsi_march:peers", "PEERS stresses not in P_k+1.");
Ft = sp0.translated_field_between(moved, new, k + 1, k + 1, F, [0; 0], extension_degree = q);
vals_t = kron(eye(4), Dk1)*Ft;
vals_t = rotate_tensor(vals_t, R);
St = sym_projection(ops, stress_moments(new, k, compliance(vals_t, c1, c2), rs, w));
end

function [e_s, e_v, n_s, n_v] = errors(mesh, k, S, V, ex, t)
% L2 errors and exact norms of sigma and v at time t
[e_s, n_s] = tensor_error(mesh, k, S, @(x, y) ex.sigma(x, y, t));
[e_v, n_v] = vector_error(mesh, k, V, @(x, y) ex.v(x, y, t));
end
